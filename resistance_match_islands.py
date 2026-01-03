from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path
from multiprocessing import Pool, cpu_count
import re
import csv

def read_fasta_sequences(faa_path):
    return list(SeqIO.parse(faa_path, "fasta"))


def compare_sequences(seqs1, seqs2, aligner, identity_threshold=0.9):
    high_matches = []
    low_matches = []
    matched_ids1 = set()
    matched_ids2 = set()

    for record1 in seqs1:
        best_match = None
        best_identity = 0

        for record2 in seqs2:
            score = aligner.score(str(record1.seq), str(record2.seq))
            max_len = max(len(record1.seq), len(record2.seq))
            identity = score / max_len if max_len > 0 else 0

            if identity > best_identity:
                best_identity = identity
                best_match = (record1.id, record2.id, identity)

        if best_match:
            if best_identity >= identity_threshold:
                high_matches.append(best_match)
                matched_ids1.add(best_match[0])
                matched_ids2.add(best_match[1])
            else:
                low_matches.append(best_match)

    unmatched1 = [r.id for r in seqs1 if r.id not in matched_ids1]
    unmatched2 = [r.id for r in seqs2 if r.id not in matched_ids2]

    return high_matches, low_matches, unmatched1, unmatched2


def process_pair(args):
    key1, key2, faa_files = args

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0

    seqs1 = read_fasta_sequences(faa_files[key1])
    seqs2 = read_fasta_sequences(faa_files[key2])

    high, low, unmatched1, unmatched2 = compare_sequences(seqs1, seqs2, aligner)

    all_identities = [m[2] for m in high + low]
    avg_identity = (sum(all_identities) / len(all_identities) * 100) if all_identities else 0.0

    summary = (
        f"{key1} vs {key2}: average similarity of {avg_identity:.2f}%.\n"
        f"- {len(unmatched1)} gene(s) from {key1} without ≥90% match in {key2}\n"
        f"- {len(unmatched2)} gene(s) from {key2} without ≥90% match in {key1}"
    )

    return (key1, key2, summary, high, low)


def extract_similarity_network(summary_file, output_csv):
    pattern = r"(.+\.faa) vs (.+\.faa): average similarity of ([\d\.]+)%"

    with open(summary_file, "r", encoding="utf-8") as infile, \
         open(output_csv, "w", newline="", encoding="utf-8") as outfile:

        writer = csv.writer(outfile)
        writer.writerow(["Source", "Target", "Weight", "Group"])

        for line in infile:
            match = re.search(pattern, line)
            if match:
                source = match.group(1).strip()
                target = match.group(2).strip()
                weight = float(match.group(3))

                source_prefix = Path(source).parts[0]
                target_prefix = Path(target).parts[0]

                if source_prefix == target_prefix:
                    continue

                if weight >= 90:
                    group = "High (>=90%)"
                elif 50 <= weight < 90:
                    group = "Medium (50-89%)"
                else:
                    group = "Low (<50%)"

                writer.writerow([source, target, weight, group])


def simplify_island_name(full_name):
    strain_match = re.search(r"(.+?)_Islands_fisher", full_name)
    island_match = re.search(r"(Resistance|Resistance)_Island_(\d+)", full_name)

    if strain_match and island_match:
        return f"{strain_match.group(1)} ({island_match.group(2)})"
    return full_name


def rename_network_nodes(input_csv):
    input_path = Path(input_csv)
    output_path = input_path.with_name("renamed_" + input_path.stem + ".csv")

    with open(input_path, newline='', encoding="utf-8") as infile, \
         open(output_path, "w", newline='', encoding="utf-8") as outfile:

        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()

        for row in reader:
            row["Source"] = simplify_island_name(row["Source"])
            row["Target"] = simplify_island_name(row["Target"])
            writer.writerow(row)

    return output_path


def generate_global_summary(base_dir, network_csv, output_txt):
    base_dir = Path(base_dir)

    islands_per_strain = {}
    total_islands = 0

    for folder in base_dir.iterdir():
        if folder.is_dir():
            aa_dir = folder / "Amino_acids"
            if aa_dir.exists():
                islands = list(aa_dir.glob("Resistance_Island_*.faa"))
                if islands:
                    islands_per_strain[folder.name] = len(islands)
                    total_islands += len(islands)

    similarity_counts = {"High": 0, "Medium": 0, "Low": 0}
    total_edges = 0

    with open(network_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_edges += 1
            if "High" in row["Group"]:
                similarity_counts["High"] += 1
            elif "Medium" in row["Group"]:
                similarity_counts["Medium"] += 1
            elif "Low" in row["Group"]:
                similarity_counts["Low"] += 1

    with open(output_txt, "w", encoding="utf-8") as out:
        out.write("GLOBAL SUMMARY OF RESISTANCE ISLANDS ANALYSIS\n")
        out.write("=" * 55 + "\n\n")

        out.write(f"Total number of strains analyzed: {len(islands_per_strain)}\n")
        out.write(f"Total number of resistance islands: {total_islands}\n\n")

        out.write("resistance islands per strain:\n")
        for strain, count in sorted(islands_per_strain.items()):
            out.write(f"- {strain}: {count} islands\n")

        out.write("\nSimilarity network summary:\n")
        out.write(f"- Total similarity edges: {total_edges}\n")
        out.write(f"- High similarity (>=90%): {similarity_counts['High']}\n")
        out.write(f"- Medium similarity (50-89%): {similarity_counts['Medium']}\n")
        out.write(f"- Low similarity (<50%): {similarity_counts['Low']}\n")


def main():
    base_dir = Path.cwd()
    faa_files = {}

    for folder in base_dir.iterdir():
        aa_dir = folder / "Amino_acids"
        if aa_dir.exists():
            for faa_file in aa_dir.glob("Resistance_Island_*.faa"):
                faa_files[f"{folder.name}/{faa_file.name}"] = faa_file

    keys = list(faa_files.keys())
    tasks = [(keys[i], keys[j], faa_files)
             for i in range(len(keys))
             for j in range(i + 1, len(keys))]

    print(f"Starting parallel processing using {cpu_count()} CPU cores")
    print(f"Total pairwise comparisons: {len(tasks)}")

    with Pool(cpu_count()) as pool:
        results = pool.map(process_pair, tasks)

    summary_file = "comparison_results_summary.txt"
    network_file = "cytoscape_network_classified.csv"
    global_summary = "global_analysis_summary.txt"

    with open(summary_file, "w", encoding="utf-8") as f:
        for _, _, summary, _, _ in results:
            f.write(summary + "\n\n")

    extract_similarity_network(summary_file, network_file)
    renamed_network = rename_network_nodes(network_file)

    generate_global_summary(
        base_dir=base_dir,
        network_csv=renamed_network,
        output_txt=global_summary
    )

    print("\nPipeline completed successfully.\nThanks for using it! \nIf you have any questions, please contact me by email at janainedpaula@gmail.com.")
    print("\nGenerated files:")
    print(f"- {summary_file}")
    print(f"- {network_file}")
    print(f"- {renamed_network}")
    print(f"- {global_summary}")


if __name__ == "__main__":
    main()
