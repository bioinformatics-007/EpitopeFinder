# plotting.py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

def plot_bcell_epitopes_from_csv(path):
    """Plot B-cell epitope prediction scores as a bar chart."""
    df = pd.read_csv(path)
    peptides = df.iloc[:, 0]
    scores = df.iloc[:, 1]
    plt.figure(figsize=(8, 4))
    plt.bar(peptides, scores, color='skyblue')
    plt.title('B-cell Epitope Prediction')
    plt.ylabel('Prediction Score')
    plt.xlabel('Peptides')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_mhci_heatmap_from_csv(path):
    """Plot MHC Class I binding affinity as a heatmap."""
    df = pd.read_csv(path, index_col=0)
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, annot=True, fmt=".1f", cmap='YlGnBu', cbar_kws={'label': 'IC50 (nM)'})
    plt.title('MHC Class I Binding Affinity (Heatmap)')
    plt.xlabel('Peptides')
    plt.ylabel('HLA Alleles')
    plt.tight_layout()
    plt.show()

def plot_mhcii_heatmap_from_csv(path):
    """Plot MHC Class II binding affinity as a heatmap."""
    df = pd.read_csv(path, index_col=0)
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, annot=True, fmt=".1f", cmap='YlOrRd', cbar_kws={'label': 'Binding Score'})
    plt.title('MHC Class II Binding Affinity (Heatmap)')
    plt.xlabel('Peptides')
    plt.ylabel('HLA Alleles')
    plt.tight_layout()
    plt.show()

def plot_ctl_epitopes_from_csv(path):
    """Plot CTL epitope prediction scores as a line plot."""
    df = pd.read_csv(path)
    peptides = df.iloc[:, 0]
    scores = df.iloc[:, 1]
    plt.figure(figsize=(8, 4))
    plt.plot(peptides, scores, marker='o', color='purple')
    plt.title('CTL Epitope Prediction Score')
    plt.ylabel('Combined Score')
    plt.xlabel('Peptides')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_proteasomal_cleavage_from_txt(path):
    """Plot proteasomal cleavage prediction scores as a stem plot."""
    with open(path) as f:
        lines = f.readlines()
    cleavage_scores = []
    peptides = []
    for line in lines:
        if line.strip() and not line.startswith("#"):
            parts = line.strip().split()
            peptides.append(parts[0])
            cleavage_scores.append(float(parts[1]))
    plt.figure(figsize=(8, 4))
    plt.stem(peptides, cleavage_scores, basefmt=" ")
    plt.title('Proteasomal Cleavage Prediction')
    plt.ylabel('Cleavage Score')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_localization_from_csv(path):
    """Plot subcellular localization distribution as a pie chart."""
    df = pd.read_csv(path)
    counts = df['Localization'].value_counts()
    labels = counts.index
    values = counts.values
    colors = ['lightgreen', 'gold', 'lightcoral', 'lightskyblue']
    plt.figure(figsize=(6, 6))
    plt.pie(values, labels=labels, autopct='%1.1f%%', colors=colors[:len(labels)])
    plt.title('Subcellular Localization Distribution')
    plt.tight_layout()
    plt.show()
