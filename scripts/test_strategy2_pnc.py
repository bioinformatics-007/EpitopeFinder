#!/usr/bin/env python3
"""
test_strategy2_pnc.py — Comprehensive Automated Diagnostic & P&C Suite for Strategy 2

Tests Strategy 2 across permutations:
1. PNC_01: Raw FASTA + All 5 tools (iapred, algpred, instability, molwt, wolfpsort) + Bacteria
2. PNC_02: Raw FASTA + All 5 tools + Virus
3. PNC_03: Raw FASTA + All 5 tools + Fungi
4. PNC_04: File Upload + Subset tools (iapred, instability, molwt) + Bacteria
5. PNC_05: UniProt P0A8T7 + All 5 tools + Bacteria
6. PNC_06: Multi-sequence FASTA + All 5 tools + Bacteria
"""

import sys
import os
import json
import time
import requests
from pathlib import Path
import pandas as pd

API_BASE_URL = os.environ.get("API_BASE_URL", "http://localhost:8000/api")
POLL_INTERVAL = 3  # seconds
MAX_WAIT_TIME = 300  # seconds

SAMPLE_FASTA = """>test_prot_seq1
MAGAASPCKANLGLALALVLLVVVHMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP
"""

MULTI_SEQ_FASTA = """>seq1_outer_mem
MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGUIAAA
>seq2_cytoplasmic
MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFP
"""

PNC_TEST_CASES = [
    {
        "id": "PNC_01_RawFASTA_AllTools_Bacteria",
        "input_type": "text",
        "input_value": SAMPLE_FASTA,
        "strategy": 2,
        "pathogen_type": "bacteria",
        "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"],
        "description": "Raw FASTA + All 5 tools + Bacteria"
    },
    {
        "id": "PNC_02_RawFASTA_AllTools_Virus",
        "input_type": "text",
        "input_value": SAMPLE_FASTA,
        "strategy": 2,
        "pathogen_type": "virus",
        "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"],
        "description": "Raw FASTA + All 5 tools + Virus"
    },
    {
        "id": "PNC_03_RawFASTA_AllTools_Fungi",
        "input_type": "text",
        "input_value": SAMPLE_FASTA,
        "strategy": 2,
        "pathogen_type": "fungi",
        "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"],
        "description": "Raw FASTA + All 5 tools + Fungi"
    },
    {
        "id": "PNC_04_FileUpload_SubsetTools_Bacteria",
        "input_type": "file",
        "input_value": SAMPLE_FASTA,
        "filename": "test_upload_s2.fasta",
        "strategy": 2,
        "pathogen_type": "bacteria",
        "selected_tools": ["iapred", "instability", "molwt"],
        "description": "File Upload + Subset tools (iapred, instability, molwt) + Bacteria"
    },
    {
        "id": "PNC_05_UniProt_P0A8T7_AllTools_Bacteria",
        "input_type": "text",
        "input_value": "P0A8T7",
        "strategy": 2,
        "pathogen_type": "bacteria",
        "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"],
        "description": "UniProt P0A8T7 + All 5 tools + Bacteria"
    },
    {
        "id": "PNC_06_MultiSeqFASTA_AllTools_Bacteria",
        "input_type": "text",
        "input_value": MULTI_SEQ_FASTA,
        "strategy": 2,
        "pathogen_type": "bacteria",
        "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"],
        "description": "Multi-Sequence FASTA + All 5 tools + Bacteria"
    }
]


def submit_job(test_case):
    """Submits job via API endpoint."""
    if test_case["input_type"] == "text":
        url = f"{API_BASE_URL}/jobs/submit"
        payload = {
            "input_value": test_case["input_value"],
            "strategy": test_case["strategy"],
            "pathogen_type": test_case["pathogen_type"],
            "selected_tools": test_case["selected_tools"]
        }
        res = requests.post(url, json=payload, timeout=10)
    else:
        url = f"{API_BASE_URL}/jobs/submit-with-file"
        files = {"file": (test_case["filename"], test_case["input_value"], "text/plain")}
        data = {
            "strategy": str(test_case["strategy"]),
            "pathogen_type": test_case["pathogen_type"],
            "selected_tools": json.dumps(test_case["selected_tools"])
        }
        res = requests.post(url, files=files, data=data, timeout=10)

    res.raise_for_status()
    return res.json()["job_id"]


def poll_job(job_id):
    """Polls job status until complete or error."""
    start_time = time.time()
    url = f"{API_BASE_URL}/jobs/{job_id}/status"
    while time.time() - start_time < MAX_WAIT_TIME:
        res = requests.get(url, timeout=5)
        if res.status_code == 200:
            data = res.json()
            status = data.get("status")
            if status == "completed":
                return True, data
            elif status == "failed":
                return False, data
        time.sleep(POLL_INTERVAL)
    return False, {"error": "Timeout waiting for job"}


def diagnose_file(job_id, output_file_info):
    """Downloads and performs deep diagnostic on an output file."""
    rel_path = output_file_info["relative_path"]
    download_url = output_file_info["download_url"]
    if not download_url.startswith("http"):
        download_url = f"{API_BASE_URL.replace('/api', '')}{download_url}"

    diag = {
        "file": rel_path,
        "filename": Path(rel_path).name,
        "size_bytes": output_file_info.get("size_bytes", 0),
        "http_status": 0,
        "valid_csv": False,
        "num_rows": 0,
        "columns": [],
        "issues": []
    }

    try:
        res = requests.get(download_url, timeout=10)
        diag["http_status"] = res.status_code
        if res.status_code != 200:
            diag["issues"].append(f"HTTP Download Error {res.status_code}")
            return diag

        content = res.text
        if not content.strip():
            diag["issues"].append("Empty File (0 content bytes)")
            return diag

        if "<html" in content.lower() or "<!doctype" in content.lower():
            diag["issues"].append("File contains HTML error page instead of data")
            return diag

        if "Error." in content:
            diag["issues"].append("File contains raw error string ('Error.')")

        if Path(rel_path).suffix.lower() == ".csv":
            from io import StringIO
            df = pd.read_csv(StringIO(content))
            diag["valid_csv"] = True
            diag["num_rows"] = len(df)
            diag["columns"] = list(df.columns)

            if any(str(col).startswith("#") for col in df.columns):
                diag["issues"].append("Header contains raw '#' comment character")

            if len(df) == 0:
                diag["issues"].append("CSV has header but 0 data rows")

    except Exception as e:
        diag["issues"].append(f"Parse Exception: {e}")

    return diag


def main():
    print("🚀 Starting Strategy 2 Automated P&C Diagnostic Suite...")
    print(f"Backend API URL: {API_BASE_URL}\n")

    overall_results = []

    for test_case in PNC_TEST_CASES:
        t_id = test_case["id"]
        desc = test_case["description"]
        print(f"--------------------------------------------------")
        print(f"▶ Running {t_id}")
        print(f"  Description: {desc}")

        case_result = {
            "id": t_id,
            "description": desc,
            "job_id": None,
            "submit_success": False,
            "job_success": False,
            "execution_time_sec": 0,
            "diagnostics": [],
            "error": None
        }

        start_time = time.time()
        try:
            job_id = submit_job(test_case)
            case_result["job_id"] = job_id
            case_result["submit_success"] = True
            print(f"  Submitted Job ID: {job_id}")

            success, status_data = poll_job(job_id)
            case_result["execution_time_sec"] = round(time.time() - start_time, 1)

            if success:
                case_result["job_success"] = True
                print(f"  [SUCCESS] Job completed in {case_result['execution_time_sec']}s")

                res_url = f"{API_BASE_URL}/jobs/{job_id}/results"
                res_resp = requests.get(res_url, timeout=10).json()
                outputs = res_resp.get("outputs", [])
                print(f"  Found {len(outputs)} output files. Performing diagnostics...")

                for out_file in outputs:
                    file_diag = diagnose_file(job_id, out_file)
                    case_result["diagnostics"].append(file_diag)
                    status_icon = "✅" if not file_diag["issues"] else "⚠️"
                    print(f"    {status_icon} {file_diag['filename']} ({file_diag['size_bytes']} bytes, {file_diag['num_rows']} rows) — Issues: {file_diag['issues'] if file_diag['issues'] else 'None'}")

            else:
                case_result["error"] = status_data.get("error", "Job failed")
                print(f"  [FAILED] Job failed: {case_result['error']}")

        except Exception as e:
            case_result["error"] = str(e)
            print(f"  [EXCEPT] Exception occurred: {e}")

        overall_results.append(case_result)

    print("\n==================================================")
    print("📊 STRATEGY 2 P&C DIAGNOSTIC SUMMARY MATRIX")
    print("==================================================")
    total_cases = len(overall_results)
    passed_cases = sum(1 for c in overall_results if c["job_success"])
    print(f"Overall Progress: {passed_cases}/{total_cases} P&C Test Cases Passed\n")

    summary_rows = []
    for c in overall_results:
        files_count = len(c["diagnostics"])
        issue_files = [d['filename'] for d in c["diagnostics"] if d['issues']]
        issue_str = ", ".join(issue_files) if issue_files else "None"
        summary_rows.append({
            "Test ID": c["id"],
            "Status": "PASSED ✅" if c["job_success"] else "FAILED ❌",
            "Job ID": c["job_id"],
            "Time (s)": c["execution_time_sec"],
            "Files": files_count,
            "File Issues": issue_str
        })

    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))

    report_file = Path("/home/amity/EpitopeFinder/EpitopeFinder/scripts/strategy2_pnc_diagnosis.json")
    with open(report_file, "w") as f:
        json.dump(overall_results, f, indent=2)
    print(f"\nSaved detailed diagnostic JSON report to: {report_file}")


if __name__ == "__main__":
    main()
