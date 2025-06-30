from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path

def read_fasta_sequences(faa_path):
    return list(SeqIO.parse(faa_path, "fasta"))

def compare_sequences(seqs1, seqs2, aligner):
    high_matches = []
    low_matches = []
    matched_ids1 = set()
    matched_ids2 = set()

    for record1 in seqs1:
        best_match = None
        best_identity = 0
        best_record2 = None

        for record2 in seqs2:
            score = aligner.score(str(record1.seq), str(record2.seq))
            max_len = max(len(record1.seq), len(record2.seq))
            identity = score / max_len if max_len > 0 else 0

            if identity > best_identity:
                best_identity = identity
                best_match = (record1.id, record2.id, identity)
                best_record2 = record2.id

        if best_match:
            if best_identity >= 0.9:
                high_matches.append(best_match)
                matched_ids1.add(best_match[0])
                matched_ids2.add(best_match[1])
            else:
                low_matches.append(best_match)

    unmatched1 = [r.id for r in seqs1 if r.id not in matched_ids1]
    unmatched2 = [r.id for r in seqs2 if r.id not in matched_ids2]

    return high_matches, low_matches, unmatched1, unmatched2

def main():
    base_dir = Path.cwd()
    subfolders = [f for f in base_dir.iterdir() if f.is_dir()]
    faa_files = {}

    for folder in subfolders:
        aa_dir = folder / "Amino_acids"
        if aa_dir.exists():
            for faa_file in aa_dir.glob("Virulence_Island_*.faa"):
                faa_files[f"{folder.name}/{faa_file.name}"] = faa_file

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0

    results_high = []
    results_low = []
    summary_all = []

    keys = list(faa_files.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            key1, key2 = keys[i], keys[j]
            seqs1 = read_fasta_sequences(faa_files[key1])
            seqs2 = read_fasta_sequences(faa_files[key2])
            high_matches, low_matches, unmatched1, unmatched2 = compare_sequences(seqs1, seqs2, aligner)

            all_identities = [m[2] for m in high_matches + low_matches]
            avg_identity = (sum(all_identities) / len(all_identities) * 100) if all_identities else 0.0

            summary = (
                f"{key1} vs {key2}: similaridade média de {avg_identity:.2f}%.\n"
                f"- {len(unmatched1)} gene(s) de {key1} sem correspondência ≥ 90% em {key2}\n"
                f"- {len(unmatched2)} gene(s) de {key2} sem correspondência ≥ 90% em {key1}"
            )
            summary_all.append(summary)

            if high_matches:
                results_high.append(f"Comparando {key1} vs {key2}:")
                for m in high_matches:
                    results_high.append(f"  {m[0]} (de {key1}) ≈ {m[1]} (de {key2}) | Identidade: {m[2]*100:.2f}%")
                results_high.append("")

            if low_matches:
                results_low.append(f"Comparando {key1} vs {key2}:")
                for m in low_matches:
                    results_low.append(f"  {m[0]} (de {key1}) ≈ {m[1]} (de {key2}) | Identidade: {m[2]*100:.2f}%")
                results_low.append("")

    with open("virulence_island_comparison_results.txt", "w") as f:
        f.write("RESUMO GERAL DAS COMPARAÇÕES ENTRE ILHAS DE VIRULÊNCIA:\n\n")
        f.write("\n\n".join(summary_all))
        f.write("\n\n===============================\n\n")
        f.write("\n".join(results_high))

    with open("low_identity_matches.txt", "w") as f:
        f.write("\n".join(results_low))

    print("Comparação finalizada.")
    print("Resultados salvos em:")
    print("- virulence_island_comparison_results.txt (resumo + matches ≥ 90%)")
    print("- low_identity_matches.txt (matches < 90%)")

main()
