import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_epitope_analysis(results_dir, uniprot_ids, epitope_types, timestamp_dir, batch_no):
    """
    Create a side-by-side bar plot in PNG format for downstream predictions across epitope types:
    - AlgPred: Allergen vs. Non-Allergen epitopes
    - IAPred: High, Moderate, and Low antigenicity epitopes
    - ToxinPred3: Toxic vs. Non-Toxic epitopes
    Each epitope type has three groups of bars (one group per tool), with AlgPred and ToxinPred3 having two bars,
    and IAPred having three bars, placed side by side.
    """
    # Initialize counts: {epitope_type: {tool: {category: count}}}
    counts = {
        et: {
            "AlgPred": {"Allergen": 0, "Non-Allergen": 0},
            "IAPred": {"High": 0, "Moderate": 0, "Low": 0},
            "ToxinPred3": {"Toxic": 0, "Non-Toxic": 0}
        } for et in epitope_types
    }

    for uniprot_id in uniprot_ids:
        for epitope_type in epitope_types:
            # File paths from strategy_5 downstream outputs
            alg_file = Path(results_dir) / f"{epitope_type}_Combined_AlgPred_{batch_no}.csv"
            iap_file = Path(results_dir) / f"{epitope_type}_combined_iapred.csv"
            tox_file = Path(results_dir) / f"{epitope_type}_Combined_ToxinPred3_{batch_no}.csv"

            # AlgPred: Count Allergen and Non-Allergen
            if alg_file.exists():
                try:
                    alg_df = pd.read_csv(alg_file, dtype=str)
                    if "Prediction" in alg_df.columns:
                        alg_df["Prediction"] = alg_df["Prediction"].str.lower()
                        counts[epitope_type]["AlgPred"]["Allergen"] += len(alg_df[alg_df["Prediction"] == "allergen"])
                        counts[epitope_type]["AlgPred"]["Non-Allergen"] += len(alg_df[alg_df["Prediction"] == "non-allergen"])
                except Exception as e:
                    print(f"AlgPred read failed for {alg_file}: {e}")

            # IAPred: Count High, Moderate, and Low
            if iap_file.exists():
                try:
                    iap_df = pd.read_csv(iap_file, dtype=str)
                    if "Antigenicity_Category" in iap_df.columns:
                        iap_df["Antigenicity_Category"] = iap_df["Antigenicity_Category"].str.lower()
                        counts[epitope_type]["IAPred"]["High"] += len(iap_df[iap_df["Antigenicity_Category"] == "high"])
                        counts[epitope_type]["IAPred"]["Moderate"] += len(iap_df[iap_df["Antigenicity_Category"] == "moderate"])
                        counts[epitope_type]["IAPred"]["Low"] += len(iap_df[iap_df["Antigenicity_Category"] == "low"])
                except Exception as e:
                    print(f"IAPred read failed for {iap_file}: {e}")

            # ToxinPred3: Count Toxic and Non-Toxic
            if tox_file.exists():
                try:
                    tox_df = pd.read_csv(tox_file, dtype=str)
                    if "Prediction" in tox_df.columns:
                        tox_df["Prediction"] = tox_df["Prediction"].str.lower()
                        counts[epitope_type]["ToxinPred3"]["Toxic"] += len(tox_df[tox_df["Prediction"].isin(["toxic", "toxin"])])
                        counts[epitope_type]["ToxinPred3"]["Non-Toxic"] += len(tox_df[tox_df["Prediction"].isin(["non-toxic", "non-toxin"])])
                except Exception as e:
                    print(f"ToxinPred3 read failed for {tox_file}: {e}")

    print(f"DEBUG: Final counts for batch {batch_no}: {counts}")

    # Prepare plot data
    tools = ["AlgPred", "IAPred", "ToxinPred3"]
    categories = {
        "AlgPred": ["Allergen", "Non-Allergen"],
        "IAPred": ["High", "Moderate", "Low"],
        "ToxinPred3": ["Toxic", "Non-Toxic"]
    }
    data = {
        tool: {
            cat: [counts[et][tool][cat] for et in epitope_types]
            for cat in categories[tool]
        } for tool in tools
    }

    total = {
        tool: sum(sum(data[tool][cat]) for cat in categories[tool])
        for tool in tools
    }
    total_count = sum(total.values())
    if total_count == 0:
        print(f"No data to plot for batch {batch_no}")
        return

    # Plot
    bar_width = 0.1  # Narrower bars to accommodate all categories side by side
    x = np.arange(len(epitope_types))
    labels = [et.replace("bcell", "B-Cell").replace("mhci", "MHC-I").replace("mhcii", "MHC-II") for et in epitope_types]
    colors = {
        "Allergen": "#99cc99",       # Light green
        "Non-Allergen": "#336633",   # Dark green
        "High": "#ff6666",           # Bright red
        "Moderate": "#cc3333",       # Medium red
        "Low": "#990000",            # Dark red
        "Toxic": "#fff68f",          # Light yellow
        "Non-Toxic": "#ffd700"       # Golden yellow
    }

    fig, ax = plt.subplots(figsize=(14, 6))  # Slightly wider figure to accommodate more bars
    offset = 0  # Track position for each bar
    max_height = 0  # Track maximum bar height for setting ylim
    for tool in tools:
        for cat in categories[tool]:
            bars = ax.bar(
                x + offset * bar_width,
                data[tool][cat],
                bar_width,
                label=f"{cat} ({tool})",
                color=colors[cat],
                edgecolor='black'
            )
            # Update max_height
            max_height = max(max_height, max(data[tool][cat]))
            # Annotate above each bar with the count
            for bar, value in zip(bars, data[tool][cat]):
                if value > 0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_y() + bar.get_height() + 0.02 * max_height,  # Reduced offset for tighter placement
                        str(value),
                        ha='center',
                        va='bottom',
                        fontsize=8
                    )
            offset += 1  # Move to next bar position

    # Increase top margin and set dynamic ylim
    plt.margins(y=0.15)  # Add 15% margin to the top
    ax.set_ylim(0, max_height * 1.2)  # Extend y-axis to 120% of max height to ensure annotations fit

    ax.set_xlabel("Epitope Type")
    ax.set_ylabel("Count of Epitopes")
    ax.set_title(f"Epitope Prediction Analysis - Batch {batch_no}")
    ax.set_xticks(x + bar_width * (offset - 1) / 2)
    ax.set_xticklabels(labels)
    ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1))
    ax.grid(True, axis='y', linestyle='--', alpha=0.5)

    plt.tight_layout()
    output_file = Path(timestamp_dir) / f"batch_{batch_no}_epitope_analysis.png"
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

    print(f"Successfully saved epitope analysis plot to {output_file}")
