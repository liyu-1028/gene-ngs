import os
import json
import logging
import pandas as pd
from dotenv import load_dotenv
from src.stage1.llm_extractor import DeepSeekExtractor

# 加载 .env 文件
load_dotenv()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_extraction_pipeline():
    MD_DIR = "data/markdown"
    PROCESSED_DIR = "data/processed"
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    
    # 替换为您的 DeepSeek API Key
    DEEPSEEK_API_KEY = os.environ.get("DEEPSEEK_API_KEY", "your_api_key_here")
    
    extractor = DeepSeekExtractor(api_key=DEEPSEEK_API_KEY)
    
    md_files = [f for f in os.listdir(MD_DIR) if f.endswith(".md")]
    if not md_files:
        logging.warning("No markdown files found to process.")
        return
        
    structured_data = []
    
    for md_filename in md_files:
        md_path = os.path.join(MD_DIR, md_filename)
        logging.info(f"Extracting features from {md_filename}...")
        
        with open(md_path, "r", encoding="utf-8") as f:
            md_content = f.read()
            
        try:
            # 限制发送给大模型的文本长度以节约 Token (根据需要调整)
            # 这里取前8000个字符，通常包含了检测报告的首页、基本信息、主要变异和小结部分
            truncated_md = md_content[:8000] 
            
            result_json = extractor.extract_features(truncated_md)
            result_json["source_file"] = md_filename
            structured_data.append(result_json)
            
            # 保存单份提取结果备查
            out_json_path = os.path.join(PROCESSED_DIR, f"{md_filename}.json")
            with open(out_json_path, "w", encoding="utf-8") as out_f:
                json.dump(result_json, out_f, ensure_ascii=False, indent=2)
                
        except Exception as e:
            logging.error(f"Failed to extract {md_filename}: {e}")

    # 合并为 CSV 供 Stage 2 训练
    if structured_data:
        df = pd.json_normalize(structured_data)
        out_csv = os.path.join(PROCESSED_DIR, "clinical_features_dataset.csv")
        df.to_csv(out_csv, index=False, encoding="utf-8")
        logging.info(f"Pipeline complete! Extracted {len(structured_data)} records to {out_csv}")

if __name__ == "__main__":
    run_extraction_pipeline()
