# modules/sasa_filter.py
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from biotite.structure import sasa
from biotite.structure.io.pdb import PDBFile

def _generate_combined_dashboard(df, output_csv):
    """
    Generates a high-end, multi-panel dashboard for SASA analysis.
    Includes a Donut Chart, Violin Distribution, and Sequence Topography.
    """
    try:
        # Set professional aesthetic
        plt.style.use('default')
        sns.set_context("paper", font_scale=1.2)
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])

        # Modern Colors
        colors = ['#2ecc71', '#e74c3c'] # Exposed Green, Buried Red
        sns_palette = sns.color_palette(["#2ecc71", "#e74c3c"])

        # 1. DONUT CHART (Top Left) - Exposure Ratio
        ax_donut = fig.add_subplot(gs[0, 0])
        counts = df['Exposure_Status'].value_counts()
        labels = counts.index
        ax_donut.pie(counts, labels=labels, autopct='%1.1f%%', startangle=140, 
                     colors=colors, pctdistance=0.85, explode=[0.05]*len(counts),
                     textprops={'fontweight': 'bold'})
        center_circle = plt.Circle((0,0), 0.70, fc='white')
        ax_donut.add_artist(center_circle)
        ax_donut.set_title("Epitope Accessibility Ratio", fontsize=15, fontweight='bold', pad=10)

        # 2. VIOLIN PLOT (Top Right) - SASA Distribution
        ax_violin = fig.add_subplot(gs[0, 1])
        sns.violinplot(data=df, y='Mean_SASA', x='Exposure_Status', palette=sns_palette, ax=ax_violin, inner="quart")
        ax_violin.axhline(y=30.0, color='black', linestyle='--', alpha=0.3)
        ax_violin.set_title("SASA Distribution Density", fontsize=15, fontweight='bold')
        ax_violin.set_ylabel("Mean SASA (Å²)")
        ax_violin.set_xlabel("Classification")

        # 3. SEQUENCE TOPOGRAPHY (Bottom - Spans 2 columns)
        ax_topo = fig.add_subplot(gs[1, :])
        
        # We use a rolling mean to show a 'surface line' for the vaccine
        df_sorted = df.sort_values('Start')
        
        # Scatter for individual epitopes
        sns.scatterplot(data=df_sorted, x='Start', y='Mean_SASA', hue='Exposure_Status', 
                        palette=sns_palette, size='Exposed_Residue_Fraction', 
                        sizes=(40, 400), alpha=0.5, ax=ax_topo, edgecolor='w')
        
        # Trend line showing vaccine surface 'peaks'
        window = max(2, len(df_sorted) // 50)
        rolling_avg = df_sorted['Mean_SASA'].rolling(window=window, center=True).mean()
        ax_topo.plot(df_sorted['Start'], rolling_avg, color='black', linewidth=2, alpha=0.6, label='Exposure Trend')

        ax_topo.axhline(y=30.0, color='red', linestyle=':', linewidth=2, label='Exposure Threshold')
        ax_topo.set_title("Vaccine Architecture Exposure Map", fontsize=16, fontweight='bold')
        ax_topo.set_xlabel("Amino Acid Position", fontsize=12)
        ax_topo.set_ylabel("Mean SASA (Å²)", fontsize=12)
        ax_topo.legend(title="Status", bbox_to_anchor=(1, 1), loc='upper left')
        
        # Add a light background "heatmap" effect
        ax_topo.fill_between(df_sorted['Start'], 0, 30, color='red', alpha=0.05)
        ax_topo.fill_between(df_sorted['Start'], 30, df_sorted['Mean_SASA'].max()+10, color='green', alpha=0.05)

        # 4. FINAL TOUCHES
        etype = os.path.basename(output_csv).split('_')[0].upper()
        fig.suptitle(f"SASA EXPOSURE DASHBOARD: {etype} EPITOPES", fontsize=22, fontweight='bold', y=0.98)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path = output_csv.replace('.csv', '_dashboard.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        return plot_path
    except Exception as e:
        print(f"Visualization error: {e}")
        return None

def run_sasa_analysis(pdb_path, epitope_csv, output_csv, sasa_threshold=30.0, fraction_threshold=0.3):
    """
    Finalized SASA analysis for epitope exposure validation with Dashboard Visualization.
    """
    if not os.path.exists(pdb_path):
        return {"status": "error", "message": f"PDB file not found: {pdb_path}"}

    try:
        # 1. Load PDB and clean structure
        pdb = PDBFile.read(pdb_path)
        structure = pdb.get_structure(model=1)
        structure = structure[structure.element != "H"]

        # 2. Calculate SASA (per atom)
        atom_sasa = sasa(structure)

        # 3. Convert to per-residue SASA
        residue_sasa = {}
        for atom, area in zip(structure, atom_sasa):
            key = (atom.chain_id, atom.res_id)
            residue_sasa.setdefault(key, 0.0)
            residue_sasa[key] += area

        # 4. Chain Selection
        available_chains = list(set([k[0] for k in residue_sasa.keys()]))
        default_chain = 'A' if 'A' in available_chains else (available_chains[0] if available_chains else 'A')

        # 5. Load Epitope CSV
        df = pd.read_csv(epitope_csv)
        if 'Start' not in df.columns or 'End' not in df.columns:
            return {"status": "error", "message": "Input CSV must contain 'Start' and 'End' columns"}

        # 6. Metrics Calculation
        def get_metrics(row):
            try:
                start, end = int(row['Start']), int(row['End'])
                vals = [residue_sasa.get((default_chain, r), 0.0) for r in range(start, end + 1)]
                if not vals: return 0.0, 0.0
                mean_val = sum(vals) / len(vals)
                exposed_count = sum(1 for v in vals if v >= sasa_threshold)
                return mean_val, exposed_count / len(vals)
            except: return 0.0, 0.0

        res_metrics = df.apply(get_metrics, axis=1)
        df['Mean_SASA'] = [m[0] for m in res_metrics]
        df['Exposed_Residue_Fraction'] = [m[1] for m in res_metrics]
        
        # 7. Classification
        df['Exposure_Status'] = df.apply(
            lambda row: 'Exposed' if (row['Mean_SASA'] >= sasa_threshold and 
                                     row['Exposed_Residue_Fraction'] >= fraction_threshold) 
            else 'Buried', axis=1
        )

        # 8. Save Data
        df.to_csv(output_csv, index=False)
        
        # 9. GENERATE COMBINED DASHBOARD
        plot_path = _generate_combined_dashboard(df, output_csv)
        
        exposed_count = len(df[df['Exposure_Status'] == 'Exposed'])
        
        return {
            "status": "success", 
            "exposed_count": exposed_count,
            "total_count": len(df),
            "chain_used": default_chain,
            "plot_saved_to": plot_path
        }

    except Exception as e:
        return {"status": "error", "message": str(e)}
