import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_file):
    df = pd.read_csv(csv_file)

    # Garantir valores numéricos
    df["Coverage_A_to_B"] = pd.to_numeric(df["Coverage_A_to_B"])
    df["Coverage_B_to_A"] = pd.to_numeric(df["Coverage_B_to_A"])

    plt.figure(figsize=(7, 6))

    plt.hexbin(
        df["Coverage_A_to_B"],
        df["Coverage_B_to_A"],
        gridsize=50,
        mincnt=1
    )

    plt.xlabel("Fraction of genes from island A covered by island B")
    plt.ylabel("Fraction of genes from island B covered by island A")
    plt.title("Directional coverage between pathogenicity islands")

    cb = plt.colorbar()
    cb.set_label("Number of island pairs")

    plt.tight_layout()
    plt.savefig("Figure_directional_coverage_hexbin.png", dpi=300)
    plt.savefig("Figure_directional_coverage_hexbin.pdf", dpi=300)

    plt.close()

    print("Figure saved as Figure_directional_coverage_hexbin.png")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: plot_directional_coverage.py <island_pairwise_coverage.csv>")
        sys.exit(1)

    main(sys.argv[1])
