import pandas as pd
import networkx as nx
import community as community_louvain
from pathlib import Path
import sys

script_dir = Path(__file__).resolve().parent

input_file = script_dir / "renamed_cytoscape_network_classified.csv"

source_col = "Source"
target_col = "Target"
weight_column = "Weight"
delimiter = ","

output_nodes = script_dir / "nodes_with_clusters.csv"
output_edges = script_dir / "edges_for_cytoscape.csv"

if not input_file.exists():
    print(f"ERROR: input file not found: {input_file}")
    sys.exit(1)

print(f"Reading input file: {input_file}")

df = pd.read_csv(input_file, sep=delimiter)
print("Detected columns:", df.columns.tolist())

G = nx.Graph()

for _, row in df.iterrows():
    G.add_edge(
        row[source_col],
        row[target_col],
        weight=row[weight_column]
    )

print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")

partition = community_louvain.best_partition(G, weight=weight_column)

nodes_df = (
    pd.DataFrame.from_dict(partition, orient="index", columns=["cluster"])
    .reset_index()
    .rename(columns={"index": "node"})
)

edges_df = df.copy()

nodes_df.to_csv(output_nodes, index=False)
edges_df.to_csv(output_edges, index=False)

print("Clustering successfully completed!")
print(f"- Nodes with cluster assignments: {output_nodes}")
print(f"- Edges for Cytoscape: {output_edges}")
