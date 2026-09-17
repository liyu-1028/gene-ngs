import os
import requests
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_single_file():
    endpoint = "http://127.0.0.1:8000/file_parse"
    # Choose a file without special characters first
    test_file = "data/raw/WS 599.4—2018.pdf"
    
    if not os.path.exists(test_file):
        logging.error(f"Test file not found: {test_file}")
        return

    payload = {
        'backend': 'pipeline',
        'formula_enable': 'true',
        'table_enable': 'true',
        'return_md': 'true',
        'return_images': 'false',
        'response_format_zip': 'false'
    }

    logging.info(f"Testing single file: {test_file}")
    with open(test_file, 'rb') as f:
        files = [('files', (os.path.basename(test_file), f, 'application/pdf'))]
        response = requests.post(endpoint, data=payload, files=files)
        
    logging.info(f"Status Code: {response.status_code}")
    try:
        logging.info(f"Response: {response.json()}")
    except:
        logging.info(f"Response Text: {response.text}")

if __name__ == "__main__":
    test_single_file()
