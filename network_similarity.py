import re
import csv
from pathlib import Path

def extract_similarity_edges(filepath, output_csv):
    # Padrão para capturar: arquivo1.faa vs arquivo2.faa: similaridade média de 95.23%
    pattern = r"(.+\.faa) vs (.+\.faa): similaridade média de ([\d\.]+)%"

    # Abre os arquivos com encoding UTF-8
    with open(filepath, "r", encoding="utf-8") as infile, open(output_csv, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Source", "Target", "Weight", "Group"])

        for line in infile:
            match = re.search(pattern, line)
            if match:
                source = match.group(1).strip()
                target = match.group(2).strip()
                weight = float(match.group(3))

                # Extrai o nome da linhagem (pasta raiz)
                source_prefix = Path(source).parts[0]
                target_prefix = Path(target).parts[0]

                # Ignora comparações dentro da mesma linhagem
                if source_prefix == target_prefix:
                    continue

                # Define o grupo de similaridade com base no peso
                if weight >= 90:
                    group = "High (≥90%)"
                elif 50 <= weight < 90:
                    group = "Medium (50–89%)"
                else:
                    group = "Low (<50%)"

                # Escreve linha no CSV
                writer.writerow([source, target, weight, group])

    print(f"Rede exportada com sucesso para '{output_csv}'.")

# Executa o script
extract_similarity_edges(
    "virulence_island_comparison_results.txt",
    "cytoscape_virulence_network_filtered.csv"
)
