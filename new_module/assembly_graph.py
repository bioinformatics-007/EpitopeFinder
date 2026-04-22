# modules/assembly_graph.py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

def plot_vaccine_architecture(n_term, components, order, linkers, c_term, output_path):
    """
    Draws a professional, high-quality schematic of the vaccine architecture.
    """
    # Modern Professional Palette
    colors = {
        "Adjuvant": "#2E86C1",  # Strong Blue
        "Linker": "#BDC3C7",     # Soft Grey
        "B-cell": "#28B463",     # Emerald Green
        "CTL": "#CB4335",        # Soft Red
        "HTL": "#D35400",        # Pumpkin Orange
        "His-tag": "#884EA0"     # Amethyst Purple
    }
    
    names = {0: "B-cell", 1: "CTL", 2: "HTL"}
    linker_names = {0: linkers[0], 1: linkers[1], 2: linkers[2]}
    
    fig, ax = plt.subplots(figsize=(14, 5))
    
    current_x = 0
    y = 0.5
    box_h = 0.3
    spacing = 0.4
    
    # helper to draw rounded boxes
    def draw_box(x, width, label, color, text_color='white', subtitle=None):
        # Rounded Rectangle
        rect = FancyBboxPatch((x, y - box_h/2), width, box_h, 
                             boxstyle="round,pad=0.02,rounding_size=0.05",
                             color=color, ec="black", lw=1)
        ax.add_patch(rect)
        
        # Primary Label
        ax.text(x + width/2, y, label, color=text_color, ha='center', va='center', 
                fontweight='bold', fontsize=10)
        
        # Subtitle (Sequence/Type)
        if subtitle:
            ax.text(x + width/2, y - 0.25, subtitle, color='black', ha='center', 
                    va='top', fontsize=8, style='italic')

    # 1. N-terminus Label
    ax.text(-0.5, y, "N-term", ha='right', va='center', fontweight='bold', color='#5D6D7E')

    # 2. Add Adjuvant
    if n_term:
        draw_box(current_x, 2.0, "Adjuvant", colors["Adjuvant"], subtitle="L7/L12 or Custom")
        current_x += 2.0 + spacing
        
        # Connection Arrow
        ax.annotate('', xy=(current_x, y), xytext=(current_x-spacing, y),
                    arrowprops=dict(arrowstyle='->', lw=1.5, color='grey'))
        
        # EAAAK Linker
        draw_box(current_x, 0.8, "EAAAK", colors["Linker"], text_color='black')
        current_x += 0.8 + spacing
        ax.annotate('', xy=(current_x, y), xytext=(current_x-spacing, y),
                    arrowprops=dict(arrowstyle='->', lw=1.5, color='grey'))

    # 3. Add Epitope Domains in selected order
    for idx in order:
        cat_name = names[idx]
        l_name = linker_names[idx]
        
        # Draw Epitope Group
        draw_box(current_x, 1.8, f"{cat_name} Epitopes", colors[cat_name], subtitle=f"Joined by {l_name}")
        current_x += 1.8 + spacing
        
        # Add arrow if not at end
        if idx != order[-1] or c_term:
            ax.annotate('', xy=(current_x, y), xytext=(current_x-spacing, y),
                        arrowprops=dict(arrowstyle='->', lw=1.5, color='grey'))

    # 4. Add His-tag
    if c_term:
        draw_box(current_x, 1.0, "His-tag", colors["His-tag"], subtitle="6x-His")
        current_x += 1.0 + spacing

    # 5. C-terminus Label
    ax.text(current_x - spacing + 0.2, y, "C-term", ha='left', va='center', fontweight='bold', color='#5D6D7E')

    # Formatting
    ax.set_xlim(-1, current_x + 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    plt.title("Constructed Multi-Epitope Vaccine Architecture", fontsize=16, fontweight='bold', pad=20)
    
    # Add a legend for clarity
    legend_patches = [mpatches.Patch(color=c, label=k) for k, c in colors.items()]
    plt.legend(handles=legend_patches, loc='lower center', bbox_to_anchor=(0.5, -0.1), 
               ncol=6, fontsize=9, frameon=False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
