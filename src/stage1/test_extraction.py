import os
import json
import logging
from dotenv import load_dotenv
from src.stage1.llm_extractor import DeepSeekExtractor

load_dotenv()
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_single_extraction():
    MD_FILE = "data/markdown/PISFZ_2310198555-1(1).md"
    
    if not os.path.exists(MD_FILE):
        logging.error(f"Markdown file not found: {MD_FILE}")
        return

    api_key = os.environ.get("DEEPSEEK_API_KEY")
    if not api_key or "your_actual_key" in api_key:
        logging.error("Valid DEEPSEEK_API_KEY not found in .env")
        return

    extractor = DeepSeekExtractor(api_key=api_key)
    
    with open(MD_FILE, "r", encoding="utf-8") as f:
        content = f.read()
        
    logging.info(f"Testing extraction on {MD_FILE} (truncated to 8000 chars)...")
    try:
        result = extractor.extract_features(content[:8000])
        print(json.dumps(result, ensure_ascii=False, indent=2))
    except Exception as e:
        logging.error(f"Extraction failed: {e}")

if __name__ == "__main__":
    test_single_extraction()
