#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clinical_auditor_v2.py  (Stage 3 v2, MSK-IMPACT 50K 量级)

相对原 src/stage3/clinical_auditor.py (root 冻结, 逻辑在此继承复用) 的变更:

1. TMB-MSI 规则升级为双向交叉校验 (原规则只有单向):
     - TMB > 20  且 MSS      -> REJECT (超高 TMB 应伴随 MSI-H/POLE, 原规则)
     - MSI-H     且 TMB < 2  -> REJECT (新增反向: MSI-H 但 TMB 极低, 逻辑矛盾)
   50K 数据源提供真实 MSI_TYPE (Stable->MSS / Instable->MSI-H)，
   无需再退化为 2017 版的 "TMB 单边校验"。
2. 两级审计以适配 5 万量级:
     Pass 1 (全量, 本地): 硬规则审计所有记录
     Pass 2 (抽样, LLM) : 随机抽取 LLM_AUDIT_N 条 (默认 2000, seed=42) 送 DeepSeek 深审
   未被抽中且通过硬规则的记录保留, 标记 audit_stage='HardRuleOnly'。
3. 输出列保持不变, 附加 audit_stage 列。

环境变量:
   LLM_AUDIT_N   LLM 深审抽样数 (默认 2000; 设为 0 可跳过 LLM)
   DEEPSEEK_API_KEY  必需 (除非 LLM_AUDIT_N=0)

用法:
   python src/stage0/clinical_auditor_v2.py \
       [--input data/synthetic/raw_synthetic_data.csv]
"""

import argparse
import ast
import json
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../stage3')))
from clinical_auditor import ClinicalAuditor   # 复用原 LLM 审计实现

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class ClinicalAuditorV2(ClinicalAuditor):
    """路由感知审计器:
      msi_mode=True  -> 双向 TMB-MSI 交叉校验 (50K 等含真实 MSI 的源)
      msi_mode=False -> TMB 单边校验 (2017 等无 MSI 字段的源, 降级模式)"""

    def __init__(self, api_key: str, msi_mode: bool = True):
        super().__init__(api_key)
        self.msi_mode = msi_mode

    def hard_rule_audit(self, row: pd.Series) -> tuple:
        tmb = row.get('biomarkers.tmb', 0)
        msi = row.get('biomarkers.msi', '') or ''
        try:
            tmb_val = float(tmb) if pd.notna(tmb) else 0.0
        except (ValueError, TypeError):
            tmb_val = 0.0

        if self.msi_mode:
            # 1a. TMB 极高但 MSS (原单向规则)
            if tmb_val > 20 and msi == 'MSS':
                return False, f"HardRule Failed: TMB extremely high ({tmb_val}) but MSI is MSS."
            # 1b. MSI-H 但 TMB 极低 (v2 新增反向规则)
            if msi == 'MSI-H' and tmb_val < 2:
                return False, f"HardRule Failed: MSI-H but TMB extremely low ({tmb_val})."
        else:
            # 单边校验降级模式: 仅合理性边界 (无 MSI 可交叉验证)
            if tmb_val < 0 or tmb_val > 400:
                return False, f"HardRule Failed: TMB out of plausible range ({tmb_val})."

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
                            return False, (f"HardRule Failed: Invalid VAF value "
                                           f"({vaf_val}) for gene {v.get('gene')}.")
        except Exception:
            pass

        return True, "Passed Hard Rules"


def run_auditor(input_csv="data/synthetic/raw_synthetic_data.csv",
                output_csv="data/synthetic/audited_synthetic_data.csv",
                log_json="data/synthetic/audit_log.json"):
    llm_n = int(os.environ.get("LLM_AUDIT_N", "2000"))
    seed = int(os.environ.get("AUDIT_SEED", "42"))

    if not os.path.exists(input_csv):
        logging.error(f"Input file not found: {input_csv}")
        return

    df = pd.read_csv(input_csv)
    logging.info(f"Loaded {len(df):,} raw synthetic records.")

    api_key = os.environ.get("DEEPSEEK_API_KEY")
    # 路由感知: msi 列全空 (如 2017 源) -> 自动降级单边校验
    msi_mode = bool((df['biomarkers.msi'].fillna('') != '').any()) if 'biomarkers.msi' in df.columns else False
    logging.info(f"MSI mode: {'cross-check (MSI available)' if msi_mode else 'single-sided TMB (no MSI in source)'}")
    auditor = ClinicalAuditorV2(api_key, msi_mode=msi_mode)

    # ---------- Pass 1: 硬规则 (全量, 本地) ----------
    logging.info("Pass 1: Hard-rule audit on ALL records ...")
    results = df.apply(lambda r: auditor.hard_rule_audit(r), axis=1)
    df['_hard_ok'] = [ok for ok, _ in results]
    hard_pass = df[df['_hard_ok']].drop(columns='_hard_ok')
    hard_reject = len(df) - len(hard_pass)
    logging.info(f"  Hard rules: PASS {len(hard_pass):,} / REJECT {hard_reject:,}")

    audit_log = [{"record_id": int(i), "status": "REJECT", "stage": "HardRule", "reason": rs}
                 for i, (ok, rs) in zip(df.index, results) if not ok]

    # ---------- Pass 2: LLM 深审 (分层随机抽样) ----------
    llm_sample_idx = pd.Index([])
    if llm_n > 0:
        if len(hard_pass) > llm_n:
            llm_sample_idx = hard_pass.sample(n=llm_n, random_state=seed).index
        else:
            llm_sample_idx = hard_pass.index
        logging.info(f"Pass 2: LLM audit on {len(llm_sample_idx):,} randomly sampled records ...")
    else:
        logging.info("Pass 2 skipped (LLM_AUDIT_N=0).")

    llm_reject_ids, llm_reasons = set(), {}

    def llm_task(index, row):
        ok, reason = auditor.llm_audit(row, index)
        return index, ok, reason

    if llm_n > 0 and llm_sample_idx.size:
        sub = hard_pass.loc[llm_sample_idx]
        done = 0
        with ThreadPoolExecutor(max_workers=5) as ex:
            futures = {ex.submit(llm_task, idx, row): idx for idx, row in sub.iterrows()}
            for fut in as_completed(futures):
                idx, ok, reason = fut.result()
                done += 1
                if not ok:
                    llm_reject_ids.add(idx)
                    llm_reasons[idx] = reason
                if done % 200 == 0:
                    logging.info(f"  LLM progress: {done:,}/{len(sub):,} "
                                 f"(rejects so far: {len(llm_reject_ids):,})")

    # ---------- 合并输出 ----------
    audited = hard_pass.drop(index=llm_reject_ids).copy()
    audited['audit_stage'] = ['LLM' if i in llm_sample_idx else 'HardRuleOnly'
                              for i in audited.index]

    for idx, row in hard_pass.iterrows():
        if idx in llm_sample_idx:
            if idx in llm_reject_ids:
                audit_log.append({"record_id": int(idx), "status": "REJECT",
                                  "stage": "LLM", "reason": llm_reasons.get(idx, ""),
                                  "data": row.to_dict()})
            else:
                audit_log.append({"record_id": int(idx), "status": "PASS",
                                  "stage": "LLM", "reason": "LLM passed"})
        else:
            audit_log.append({"record_id": int(idx), "status": "PASS",
                              "stage": "HardRuleOnly", "reason": "Not sampled for LLM"})

    audited.to_csv(output_csv, index=False, encoding="utf-8")
    with open(log_json, "w", encoding="utf-8") as f:
        json.dump(audit_log, f, ensure_ascii=False)

    logging.info(f"Audit complete. Final passed: {len(audited):,}/{len(df):,} "
                 f"(HardRule reject {hard_reject:,}, LLM reject {len(llm_reject_ids):,})")
    logging.info(f"  -> {output_csv}")
    logging.info(f"  audit_stage 分布: {audited['audit_stage'].value_counts().to_dict()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default="data/synthetic/raw_synthetic_data.csv")
    ap.add_argument('--output-csv', default="data/synthetic/audited_synthetic_data.csv")
    ap.add_argument('--log-json', default="data/synthetic/audit_log.json")
    args = ap.parse_args()
    run_auditor(input_csv=args.input, output_csv=args.output_csv, log_json=args.log_json)
