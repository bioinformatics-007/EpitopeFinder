import time
import requests

BASE_URL = "http://localhost:8000/api/jobs"
FASTA = ">sp|TEST_ID|TEST\nMANIINLWNGIVPMVQDVNVASITAFKSMIDETWDKKIEANTCISRKHRNIIHEVIRDFMKAYPKMDENRKSPLGAPMQWLTQYYILKNEYHKTMLAYDDGSLNTKFKTLNIYMITNVGQ"

strategies = [1, 2, 3, 4, 5, 6]

def submit_job(strategy):
    payload = {
        "input_value": FASTA,
        "strategy": strategy,
        "pathogen_type": "virus",
        "mhci_method": "f",
        "mhcii_method": "nmel"
    }
    if strategy == 1:
        payload["selected_tools"] = ["mhc_i", "mhc_ii", "bcell"]
    elif strategy == 3:
        payload["selected_tools"] = ["algpred2", "toxinpred2"]
        
    res = requests.post(f"{BASE_URL}/submit", json=payload)
    if res.status_code == 200:
        return res.json()["job_id"]
    else:
        print(f"Failed to submit strategy {strategy}: {res.text}")
        return None

def main():
    print("Starting deep testing of all strategies...")
    job_ids = {}
    for s in strategies:
        jid = submit_job(s)
        if jid:
            print(f"Submitted Strategy {s} -> {jid}")
            job_ids[s] = jid
    
    print("Polling jobs...")
    while job_ids:
        for s in list(job_ids.keys()):
            jid = job_ids[s]
            res = requests.get(f"{BASE_URL}/{jid}/status")
            if res.status_code == 200:
                data = res.json()
                status = data["status"]
                if status in ["completed", "failed"]:
                    print(f"Strategy {s} finished with status: {status}")
                    if status == "failed" or data.get("failed_tools"):
                        print(f"ERROR DETAILS: {data}")
                    else:
                        print(f"Strategy {s} SUCCESS: {data}")
                    del job_ids[s]
        time.sleep(3)
        
    print("All tests completed.")

if __name__ == "__main__":
    main()
