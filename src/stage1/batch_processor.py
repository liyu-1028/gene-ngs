import os
import logging
from src.stage1.mineru_client import MinerUClient

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_batch_governance():
    # 路径配置
    RAW_DIR = "data/raw"
    MD_DIR = "data/markdown"
    os.makedirs(MD_DIR, exist_ok=True)

    client = MinerUClient(endpoint_url="http://127.0.0.1:8000/file_parse")

    # 获取所有待处理的 PDF
    pdf_files = [os.path.join(RAW_DIR, f) for f in os.listdir(RAW_DIR) if f.lower().endswith(".pdf")]
    
    if not pdf_files:
        logging.warning(f"No PDF files found in {RAW_DIR}")
        return

    logging.info(f"Found {len(pdf_files)} PDF files. Starting batch processing...")

    # 每 1 个文件一批 (由于 MinerU 处理较慢，改为 1 以避免超时)
    batch_size = 1
    for i in range(0, len(pdf_files), batch_size):
        batch = pdf_files[i:i + batch_size]
        try:
            results = client.parse_pdfs_batch(batch)
            
            # 保存结果
            for file_name, md_content in results.items():
                output_path = os.path.join(MD_DIR, f"{file_name}.md")
                with open(output_path, "w", encoding="utf-8") as f:
                    f.write(md_content)
                logging.info(f"Saved: {output_path}")
                
        except Exception as e:
            logging.error(f"Failed to process batch starting with {batch[0]}: {e}")

if __name__ == "__main__":
    run_batch_governance()
