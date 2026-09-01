import requests
import time
import sys
import os

# Configuration
API_URL = "http://localhost:8000/api"
TEST_FASTA = ">test_seq\nMAGAASPCKANLGLALALVLLVVV"
MAX_POLL_RETRIES = 30
POLL_INTERVAL = 10

def run_smoke_test():
    print("--- Starting EpitopePred Smoke Test ---")
    
    # 1. Check health endpoint
    try:
        resp = requests.get(f"{API_URL}/health", timeout=5)
        if resp.status_code == 200:
            print("[OK] Backend health check passed.")
        else:
            print(f"[FAIL] Backend health check failed: {resp.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] Could not connect to API: {e}")
        return False

    # 2. Submit a job
    print("Submitting test job...")
    files = {'file': ('test.fasta', TEST_FASTA)}
    data = {
        'strategy': 1,
        'pathogen_type': 'virus',
        'mhci_method': 'f',
        'mhcii_method': 'nmel',
        'selected_tools': 'mhc1,mhc2'
    }
    
    try:
        resp = requests.post(f"{API_URL}/jobs/submit-with-file", files=files, data=data, timeout=10)
        if resp.status_code != 200:
            print(f"[FAIL] Job submission failed: {resp.text}")
            return False
        
        job_data = resp.json()
        job_id = job_data.get("job_id")
        print(f"[OK] Job submitted successfully. ID: {job_id}")
    except Exception as e:
        print(f"[FAIL] Error during submission: {e}")
        return False

    # 3. Poll for status
    print(f"Polling for job status (timeout: {MAX_POLL_RETRIES * POLL_INTERVAL}s)...")
    for i in range(MAX_POLL_RETRIES):
        try:
            resp = requests.get(f"{API_URL}/jobs/{job_id}/status", timeout=5)
            if resp.status_code != 200:
                print(f"[FAIL] Status check failed: {resp.text}")
                return False
            
            status_data = resp.json()
            status = status_data.get("status")
            print(f"  Attempt {i+1}: Status is {status}")
            
            if status == "completed":
                print("[SUCCESS] Job completed successfully!")
                return True
            if status == "failed":
                print(f"[FAIL] Job failed. Error: {status_data.get('error')}")
                return False
            
            time.sleep(POLL_INTERVAL)
        except Exception as e:
            print(f"[FAIL] Error during polling: {e}")
            return False
            
    print("[FAIL] Job timed out.")
    return False

if __name__ == "__main__":
    success = run_smoke_test()
    sys.exit(0 if success else 1)
