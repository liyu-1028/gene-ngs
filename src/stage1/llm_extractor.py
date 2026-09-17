# --- src/stage1/llm_extractor.py ---
import json
import os
from openai import OpenAI

class DeepSeekExtractor:
    def __init__(self, api_key: str, config_path: str = "configs/stage1_prompt.json"):
        self.client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        
        # 加载外部 Prompt 配置文件
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
            
        with open(config_path, "r", encoding="utf-8") as f:
            config = json.load(f)
            
        self.system_prompt = config.get("system_prompt", "")
        self.schema = config.get("extraction_schema", {})

    def extract_features(self, md_content: str) -> dict:
        response = self.client.chat.completions.create(
            model="deepseek-chat",
            messages=[
                {"role": "system", "content": self.system_prompt},
                {"role": "user", "content": f"Schema: {json.dumps(self.schema, ensure_ascii=False)}\n\nReport:\n{md_content}"}
            ],
            response_format={"type": "json_object"}
        )
        
        return json.loads(response.choices[0].message.content)