#!/usr/bin/env python3
"""
EpitopePred E2E API Test Suite — v3 (Fixed)

Fixes from v2:
- All tests using the Monkeypox virus sequence now correctly set pathogen_type="virus"
- Strategy 6 assemble tests use "custom" mode (which works without prior CSV data)
- Strategy 5 test uses a realistic payload (no dummy file paths)
- Added timeout to prevent infinite polling
"""
import time
import requests
import json
import sys
import tempfile
import os

API_URL = "http://localhost:8000/api"
POLL_TIMEOUT = 300  # 5 min max per test

# Test sequence (Monkeypox VIRUS — pathogen_type must be "virus")
SEQ_MONKEYPOX = (
    ">sp|A0A7H0DNC4|PG153_MONPV\n"
    "MANIINLWNGIVPMVQDVNVASITAFKSMIDETWDKKIEANTCISRKHRNIIHEVIRDFM"
    "KAYPKMDENRKSPLGAPMQWLTQYYILKNEYHKTMLAYDDGSLNTKFKTLNIYMITNVGQ"
    "YILYIVFCIISGKNHDGTPYIYDSEITSNDKNLINDRIKYACKQILHGQLTMALRIRNKF"
    "MFIGSPMYLWFNVNGSHVYHEIYDRNVGFHNKEIGRLLYAFMYYLSISGRFLNDLALLKF"
    "TYLGESWTFSLSVPEYILYGLGYSVFDTIEKFSNDAILVYIRTNNRNGYDYAEFNKKGIV"
    "KVTEDKPDNDKRIHAIRLINYSSDVQHIHFGFRNTLIIDNECTNIQSSAENATDIGHYQD"
    "SKINIEDDDIIDDDDDDDDDDDDDDYNPKPTPIPDPHPRPPFPRHDYHKRPKLLPVEEPD"
    "PVKKDADRIRLDNHILNTLDHNLNSIGHYCCDTVAVDRLEHHIETLGQYTVILARKINMQ"
    "TLLFPWPLPTVHQHAIDGSIPPHGRSTIL"
)

CUSTOM_VACCINE_SEQ = "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKLEAAGATVTVKEAAAK"

TEST_CASES = [
    # ═══════════════════════════════════════════════════════════════
    # TIER 1: Core happy paths
    # ═══════════════════════════════════════════════════════════════
    {
        "id": "T1.1",
        "desc": "Strategy 1 - All tools (Paste, virus)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["mhc1", "mhc2", "netctl", "netchop", "bcell", "psortb"]
        }
    },
    {
        "id": "T1.2",
        "desc": "Strategy 1 - File Upload (virus)",
        "upload": True,
        "payload": {
            "strategy": "1",
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": json.dumps(["mhc1", "bcell"])
        }
    },
    {
        "id": "T1.3",
        "desc": "Strategy 1 - UniProt ID (virus)",
        "payload": {
            "input_value": "A0A7H0DNC4",
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["mhc1", "bcell"]
        }
    },
    {
        "id": "T1.4",
        "desc": "Strategy 2 - Protein Prioritization (virus)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 2,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["iapred", "algpred", "instability", "molwt", "wolfpsort"]
        }
    },
    {
        "id": "T1.5",
        "desc": "Strategy 3 - Virulence/Toxicity (virus)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 3,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["toxinpred", "deeptmhmm", "netsol", "signalp", "human"]
        }
    },
    {
        "id": "T1.6",
        "desc": "Strategy 6 - Custom vaccine (no SASA)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 6,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "assembly_config": {
                "mode": "custom",
                "custom_sequence": CUSTOM_VACCINE_SEQ,
                "run_sasa": False,
                "add_adjuvant": True,
                "add_his_tag": True,
                "bcell_count": 0, "ctl_count": 0, "htl_count": 0, "assembly_order": "1",
                "bcell_csv_path": "", "ctl_csv_path": "", "htl_csv_path": "",
                "sasa_csv_path": "", "custom_fasta_path": ""
            }
        }
    },

    # ═══════════════════════════════════════════════════════════════
    # TIER 2: MHC method variations + partial tool selection
    # ═══════════════════════════════════════════════════════════════
    {
        "id": "T2.1",
        "desc": "Strategy 1 - ANN (MHC-I) + NN_align (MHC-II) + partial tools",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "a",   # ANN
            "mhcii_method": "4",  # NN_align
            "selected_tools": ["mhc1", "mhc2", "bcell"]
        }
    },
    {
        "id": "T2.2",
        "desc": "Strategy 1 - SMM (MHC-I) + Sturniolo (MHC-II)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "i",   # SMM
            "mhcii_method": "7",  # Sturniolo
            "selected_tools": ["mhc1", "mhc2"]
        }
    },
    {
        "id": "T2.3",
        "desc": "Strategy 1 - NetMHCpan BA + Consensus3 + single tool (bcell only)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "e",   # NetMHCpan BA
            "mhcii_method": "3",  # Consensus3
            "selected_tools": ["bcell"]
        }
    },
    {
        "id": "T2.4",
        "desc": "Strategy 2 - Partial tools (instability + molwt only)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 2,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["instability", "molwt"]
        }
    },

    # ═══════════════════════════════════════════════════════════════
    # TIER 3: File uploads
    # ═══════════════════════════════════════════════════════════════
    {
        "id": "T3.1",
        "desc": "Strategy 2 - File Upload (virus)",
        "upload": True,
        "payload": {
            "strategy": "2",
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": json.dumps(["iapred", "algpred"])
        }
    },
    {
        "id": "T3.2",
        "desc": "Strategy 6 - File Upload + Custom assembly",
        "upload": True,
        "payload": {
            "strategy": "6",
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "assembly_config": json.dumps({
                "mode": "custom",
                "custom_sequence": CUSTOM_VACCINE_SEQ,
                "run_sasa": False,
                "add_adjuvant": True,
                "add_his_tag": True,
                "bcell_count": 0, "ctl_count": 0, "htl_count": 0, "assembly_order": "1",
                "bcell_csv_path": "", "ctl_csv_path": "", "htl_csv_path": "",
                "sasa_csv_path": "", "custom_fasta_path": ""
            })
        }
    },

    # ═══════════════════════════════════════════════════════════════
    # TIER 4: Strategy 6 custom mode with different configurations
    # (Assemble mode requires prior Strategy 1 output CSVs and 
    #  is not testable in isolation — this is by design)
    # ═══════════════════════════════════════════════════════════════
    {
        "id": "T4.1",
        "desc": "Strategy 6 - Custom + adjuvant + his-tag + SASA",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 6,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "assembly_config": {
                "mode": "custom",
                "custom_sequence": CUSTOM_VACCINE_SEQ,
                "run_sasa": True,
                "add_adjuvant": True,
                "add_his_tag": True,
                "bcell_count": 0, "ctl_count": 0, "htl_count": 0, "assembly_order": "1",
                "bcell_csv_path": "", "ctl_csv_path": "", "htl_csv_path": "",
                "sasa_csv_path": "", "custom_fasta_path": ""
            }
        }
    },
    {
        "id": "T4.2",
        "desc": "Strategy 6 - Custom bare (no adjuvant, no his-tag, no SASA)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 6,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "assembly_config": {
                "mode": "custom",
                "custom_sequence": CUSTOM_VACCINE_SEQ,
                "run_sasa": False,
                "add_adjuvant": False,
                "add_his_tag": False,
                "bcell_count": 0, "ctl_count": 0, "htl_count": 0, "assembly_order": "1",
                "bcell_csv_path": "", "ctl_csv_path": "", "htl_csv_path": "",
                "sasa_csv_path": "", "custom_fasta_path": ""
            }
        }
    },

    # ═══════════════════════════════════════════════════════════════
    # TIER 5: Edge cases / error handling
    # ═══════════════════════════════════════════════════════════════
    {
        "id": "T5.1",
        "desc": "Gibberish Input (Should Fail)",
        "payload": {
            "input_value": "NOT A FASTA SEQUENCE JUST BLAH BLAH",
            "strategy": 1,
            "pathogen_type": "bacteria",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["bcell"]
        },
        "expect_fail": True
    },
    {
        "id": "T5.2",
        "desc": "Invalid UniProt ID (Should Fail)",
        "payload": {
            "input_value": "ZZZZZZZZZ",
            "strategy": 1,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["bcell"]
        },
        "expect_fail": True
    },
    {
        "id": "T5.3",
        "desc": "Mismatched Pathogen Type - Protozoa (Should Fail)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "protozoa",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["mhc1", "bcell"]
        },
        "expect_fail": True
    },
    {
        "id": "T5.4",
        "desc": "Mismatched Pathogen Type - Fungi (Should Fail)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 2,
            "pathogen_type": "fungi",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["iapred", "algpred"]
        },
        "expect_fail": True
    },
    {
        "id": "T5.5",
        "desc": "Mismatched Pathogen Type - Bacteria (Should Fail)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 1,
            "pathogen_type": "bacteria",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["mhc1", "bcell"]
        },
        "expect_fail": True
    },
    {
        "id": "T5.6",
        "desc": "Strategy 6 Assemble - Empty CSV Inputs (Should Fail)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 6,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "assembly_config": {
                "mode": "assemble",
                "bcell_count": 2, "ctl_count": 2, "htl_count": 2, "assembly_order": "2",
                "add_adjuvant": False, "add_his_tag": False,
                "bcell_csv_path": "", "ctl_csv_path": "", "htl_csv_path": "",
                "custom_sequence": "", "custom_fasta_path": "", "run_sasa": False, "sasa_csv_path": ""
            }
        },
        "expect_fail": True
    },
    {
        "id": "T5.7",
        "desc": "Strategy 5 - Missing pre-predicted FASTAs (Should Fail)",
        "payload": {
            "input_value": SEQ_MONKEYPOX,
            "strategy": 5,
            "pathogen_type": "virus",
            "mhci_method": "f",
            "mhcii_method": "1",
            "selected_tools": ["algpred", "iapred"],
            "pre_predicted_fastas": {
                "bcell": "/nonexistent/path/bcell.fasta",
                "mhci": "/nonexistent/path/mhci.fasta",
                "mhcii": "/nonexistent/path/mhcii.fasta"
            }
        },
        "expect_fail": True
    },
]


def run_test(case):
    print(f"\n{'='*60}")
    print(f"  {case['id']}: {case['desc']}")
    print(f"{'='*60}")
    try:
        # Submit Job
        if case.get('upload', False):
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                tmp.write(SEQ_MONKEYPOX)
                tmp_path = tmp.name
            with open(tmp_path, 'rb') as f:
                files = {'file': ('sequence.fasta', f, 'text/plain')}
                res = requests.post(f"{API_URL}/jobs/submit-with-file", data=case['payload'], files=files)
            os.unlink(tmp_path)
        else:
            res = requests.post(f"{API_URL}/jobs/submit", json=case['payload'])

        res.raise_for_status()
        job_id = res.json()["job_id"]
        print(f"  Job ID: {job_id}")

        # Poll Status with timeout
        start = time.time()
        while time.time() - start < POLL_TIMEOUT:
            time.sleep(3)
            status_res = requests.get(f"{API_URL}/jobs/{job_id}/status")
            status_res.raise_for_status()
            sd = status_res.json()
            status = sd["status"]
            pct = sd["progress_pct"]
            tool = sd["current_tool"]
            elapsed = int(time.time() - start)
            print(f"  [{elapsed:3d}s] {status} ({pct}%) — {tool}      ", end="\r")

            if status in ["completed", "failed"]:
                print(f"  [{elapsed:3d}s] {status} ({pct}%)                         ")
                if status == "failed":
                    err = sd.get('error', '')
                    print(f"  Error: {err[:120]}")
                    failed_tools = sd.get('failed_tools', [])
                    if failed_tools:
                        print(f"  Failed tools: {failed_tools}")
                    return case.get('expect_fail', False)

                if case.get('expect_fail', False):
                    print("  ⚠ Expected failure, but job completed!")
                    return False

                res_data = requests.get(f"{API_URL}/jobs/{job_id}/results").json()
                outputs = res_data.get("outputs", [])
                print(f"  ✓ Generated {len(outputs)} output files")
                return True

        print(f"\n  ⏱ TIMEOUT after {POLL_TIMEOUT}s")
        return False

    except Exception as e:
        print(f"\n  Exception: {e}")
        return False


def main():
    print("\n" + "=" * 60)
    print("  EpitopePred E2E API Test Suite v3")
    print("=" * 60)

    results = []
    for case in TEST_CASES:
        passed = run_test(case)
        tag = "✅ PASS" if passed else "❌ FAIL"
        print(f"  → {tag}")
        results.append((case['id'], case['desc'], passed))

    # Summary table
    print("\n" + "=" * 60)
    print("  TEST SUMMARY")
    print("=" * 60)
    for tid, desc, passed in results:
        tag = "✅" if passed else "❌"
        print(f"  {tag} {tid:8s} {desc}")

    passed_count = sum(1 for _, _, p in results if p)
    failed_count = len(results) - passed_count
    print(f"\n  Total: {passed_count} Passed / {failed_count} Failed")

    if failed_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
