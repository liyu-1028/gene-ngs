# --- src/stage1/mineru_client.py ---
import requests
import logging
import os

class MinerUClient:
    def __init__(self, endpoint_url: str = "http://127.0.0.1:8000/file_parse"):
        self.endpoint = endpoint_url

    def parse_pdfs_batch(self, pdf_file_paths: list[str]) -> dict[str, str]:
        """
        根据 MinerU 3.0.9 接口规范，批量上传 PDF (建议每批5个) 并返回文件名与 Markdown 内容的映射
        """
        logging.info(f"Sending batch of {len(pdf_file_paths)} files to MinerU at {self.endpoint}...")

        payload = {
            'backend': 'pipeline',
            'formula_enable': 'true',
            'table_enable': 'true',
            'return_md': 'true',
            'return_images': 'false',
            'response_format_zip': 'false'
        }

        files = []
        file_handles = []
        try:
            # 批量打开文件
            for path in pdf_file_paths:
                f = open(path, 'rb')
                file_handles.append(f)
                files.append(('files', (os.path.basename(path), f, 'application/pdf')))

            response = requests.post(self.endpoint, data=payload, files=files)
            response.raise_for_status()

            resp_data = response.json()
            results_map = {}

            if resp_data.get("status") == "completed":
                results = resp_data.get("results", {})
                for file_id, content in results.items():
                    if "md_content" in content:
                        results_map[file_id] = content["md_content"]

                logging.info(f"Successfully processed {len(results_map)} files in this batch.")
                return results_map
            else:
                error_msg = resp_data.get("error") or "Unknown error"
                raise RuntimeError(f"MinerU batch processing failed: {error_msg}")

        except Exception as e:
            logging.error(f"Error during MinerU batch parsing: {e}")
            raise
        finally:
            # 确保所有文件句柄关闭
            for fh in file_handles:
                fh.close()