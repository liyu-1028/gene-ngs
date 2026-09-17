import os
import json
import logging
import pandas as pd
import ast
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv
from openai import OpenAI
import time

load_dotenv()
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class ClinicalAuditor:
    def __init__(self, api_key: str):
        # 使用 deepseek-chat 来替代 deepseek-reasoner 以获得更稳定的 JSON 输出和更快的速度
        # 如果必须使用 R1，可以改回 deepseek-reasoner，但需要手动解析 JSON
        self.client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        self.audit_logs = []

    def hard_rule_audit(self, row: pd.Series) -> tuple[bool, str]:
        """硬性规则审计"""
        # 1. TMB 与 MSI 关联性检查
        tmb = row.get('biomarkers.tmb', 0)
        msi = row.get('biomarkers.msi', 'MSS')
        
        try:
            tmb_val = float(tmb) if pd.notna(tmb) else 0.0
        except ValueError:
            tmb_val = 0.0
            
        if tmb_val > 20 and msi == 'MSS':
            return False, f"HardRule Failed: TMB is extremely high ({tmb_val}) but MSI is MSS."

        # 2. VAF 丰度检查
        variants_str = row.get('somatic_variants', '[]')
        try:
            variants = ast.literal_eval(variants_str) if isinstance(variants_str, str) else variants_str
            if isinstance(variants, list):
                for v in variants:
                    vaf = v.get('vaf')
                    if vaf is not None:
                        vaf_val = float(vaf)
                        if vaf_val < 0 or vaf_val > 1.0:
                            return False, f"HardRule Failed: Invalid VAF value ({vaf_val}) for gene {v.get('gene')}."
        except Exception:
            pass # 忽略解析错误，交由后续处理
            
        return True, "Passed Hard Rules"

    def llm_audit(self, row: pd.Series, record_id: int) -> tuple[bool, str]:
        """LLM 深度逻辑审计"""
        patient_desc = (
            f"患者年龄: {row.get('patient_info.age')}, "
            f"性别: {row.get('patient_info.gender')}, "
            f"病理诊断: {row.get('patient_info.pathological_diagnosis')}, "
            f"TMB: {row.get('biomarkers.tmb')}, "
            f"MSI: {row.get('biomarkers.msi')}, "
            f"体细胞变异: {row.get('somatic_variants')}"
        )

        prompt = f"""
你是一位资深的肿瘤学与临床生物信息学专家。请评估以下虚拟生成的 NGS 测序结果是否在生物学和临床上具有合理性。
如果存在严重的医学常识错误（例如：性别与癌种矛盾、极度罕见的驱动基因共突变组合等），请拒绝。
如果基本合理（即使是罕见情况，但在理论上可能存在），请通过。

病例信息:
{patient_desc}

请严格输出 JSON 格式：
{{
    "status": "PASS" 或者 "REJECT",
    "reason": "你的详细医学推理理由"
}}
"""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                response = self.client.chat.completions.create(
                    model="deepseek-chat",
                    messages=[
                        {"role": "system", "content": "You are a helpful medical data auditor. Output valid JSON only."},
                        {"role": "user", "content": prompt}
                    ],
                    response_format={"type": "json_object"},
                    timeout=30
                )
                result_str = response.choices[0].message.content
                result = json.loads(result_str)
                is_pass = (result.get("status", "REJECT") == "PASS")
                return is_pass, result.get("reason", "No reason provided")
            except Exception as e:
                if attempt == max_retries - 1:
                    return False, f"LLM Error: {e}"
                time.sleep(2)

def run_auditor():
    INPUT_CSV = "data/synthetic/raw_synthetic_data.csv"
    OUTPUT_CSV = "data/synthetic/audited_synthetic_data.csv"
    LOG_JSON = "data/synthetic/audit_log.json"
    
    if not os.path.exists(INPUT_CSV):
        logging.error(f"Input file not found: {INPUT_CSV}")
        return

    df = pd.read_csv(INPUT_CSV)
    api_key = os.environ.get("DEEPSEEK_API_KEY")
    auditor = ClinicalAuditor(api_key=api_key)
    
    audit_logs = []
    audited_records = []
    
    sample_df = df # 全量审计
    logging.info(f"Starting audit for {len(sample_df)} records...")

    def process_record(index, row):
        record_id = index
        
        # 1. Hard Rule
        passed_hard, hr_reason = auditor.hard_rule_audit(row)
        if not passed_hard:
            return {
                "record_id": record_id,
                "status": "REJECT",
                "stage": "HardRule",
                "reason": hr_reason,
                "data": row.to_dict()
            }, None

        # 2. LLM Audit
        passed_llm, llm_reason = auditor.llm_audit(row, record_id)
        if not passed_llm:
            return {
                "record_id": record_id,
                "status": "REJECT",
                "stage": "LLM",
                "reason": llm_reason,
                "data": row.to_dict()
            }, None
            
        # Passed
        return {
            "record_id": record_id,
            "status": "PASS",
            "stage": "Final",
            "reason": llm_reason,
            "data": row.to_dict()
        }, row

    # 并发执行以加速 API 调用
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(process_record, idx, row): idx for idx, row in sample_df.iterrows()}
        
        for future in as_completed(futures):
            log_entry, passed_row = future.result()
            audit_logs.append(log_entry)
            if passed_row is not None:
                audited_records.append(passed_row)

    # 保存日志
    with open(LOG_JSON, "w", encoding="utf-8") as f:
        json.dump(audit_logs, f, ensure_ascii=False, indent=2)
        
    # 保存过滤后的数据
    if audited_records:
        audited_df = pd.DataFrame(audited_records)
        audited_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
        logging.info(f"Audit completed. Passed: {len(audited_df)}/{len(sample_df)}. Saved to {OUTPUT_CSV}")
    else:
        logging.warning("All records were rejected by the auditor.")
        
    logging.info(f"Detailed audit logs saved to {LOG_JSON}")

if __name__ == "__main__":
    run_auditor()
