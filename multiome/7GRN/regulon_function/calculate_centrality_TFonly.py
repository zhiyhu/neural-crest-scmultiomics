##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate Centrality of TF-gene networks in global GRN
## Zhiyuan Hu, ChatGPT
## 26 Apr 2023
## last modified 15 Jan 2024
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("/home/huzhiy/miniforge3/envs/scenicplus")
# Refer to https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html

import pandas as pd
import networkx as nx
import os
import matplotlib.pyplot as plt

# Read the CSV file
work_dir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/'
os.chdir(work_dir)
df = pd.read_csv('scenicplus/cytoscape/eRegulon_metadata_filtered.csv',index_col=0)
df

tf_metadata = pd.read_csv('../../data/motif_with_newm/dr11_motif2tf_curated_filtered20230122.tbl', sep='\t', index_col = 0)
# Get the list of gene names from tf_metadata
tf_names = tf_metadata['gene_name'].tolist()
# Subset df by the gene names in gene_names
df = df[df['Gene'].isin(tf_names)]

# Remove duplicates
df_unique = df.drop_duplicates(subset=['TF', 'Gene'])
df_unique

# Create a directed graph from the DataFrame
G = nx.from_pandas_edgelist(df, source="TF", target="Gene", edge_attr=["TF2G_regulation"], create_using=nx.DiGraph())

#~~~~~~~~~~~~~~~~~~#
# Modify the Graph #
#~~~~~~~~~~~~~~~~~~#

# Define a function to map values to colors
def map_to_color(value):
    if value == "activation":
        return "green"
    elif value == "repression":
        return "red"
    else:
        return "black"

# Modify the arrow direction and color based on the "TF2G_regulation" column
edge_colors = [map_to_color(value) for value in nx.get_edge_attributes(G, "TF2G_regulation").values()]
edge_arrows = ["from" if value == "activation" else "to" for value in nx.get_edge_attributes(G, "TF2G_regulation").values()]

# Draw the graph
plt.figure(figsize=(12,12))
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_color="lightblue")
nx.draw_networkx_edges(G, pos, arrowsize=20, edge_color=edge_colors, arrowstyle="-|>", connectionstyle="arc3,rad=0.1", arrows=edge_arrows)
nx.draw_networkx_labels(G, pos)
plt.axis("off")
plt.show()

os.chdir('multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/figures/subgraph')
tf_selected = ['fli1a','rxraa','rarga','erfl3','elk3','nr2f5']
for tf_iter in tf_selected:
    # Subset the graph by nodes 'mite' and 'tfec' and their neighbors within 1 step
    subgraph = nx.ego_graph(G,  tf_iter, radius=1, center=True, undirected=True)

    # Get the list of mutual edges
    mutual_edges = [(u, v) for u, v, d in subgraph.edges(data=True) if (v, u) in subgraph.edges() and u != v]

    # Draw the graph
    plt.figure(figsize=(12,12))
    pos = nx.circular_layout(subgraph)
    nx.draw_networkx_nodes(subgraph, pos, node_color="lightblue")
    nx.draw_networkx_edges(subgraph, pos, arrowstyle='->',  arrowsize=10, width=1, alpha=0.5, arrows=edge_arrows)
    nx.draw_networkx_edges(subgraph, pos, edgelist=mutual_edges, arrowstyle='->', arrowsize=10, width=1, alpha=0.5, edge_color='red')
    nx.draw_networkx_labels(subgraph, pos, font_size=8)
    plt.axis("off")
    plt.savefig(f'{tf_iter}_subset_graph.png', dpi=300, bbox_inches='tight')
    plt.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compute different centrality measures #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
in_degree_centrality = nx.in_degree_centrality(G)
out_degree_centrality = nx.out_degree_centrality(G)
closeness_centrality = nx.closeness_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)
eigenvector_centrality = nx.eigenvector_centrality_numpy(G)

def top_n_centrality(centrality, n=10):
    # Sort the centrality dictionary by values (largest to smallest)
    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)

    # Return the top n entries
    return sorted_centrality[:n]

# Print the top 10 nodes for each centrality measure
print("Top 10 In-degree centrality:", top_n_centrality(in_degree_centrality))
print("Top 10 Out-degree centrality:", top_n_centrality(out_degree_centrality))
print("Top 10 Closeness centrality:", top_n_centrality(closeness_centrality))
print("Top 10 Betweenness centrality:", top_n_centrality(betweenness_centrality))
print("Top 10 Eigenvector centrality:", top_n_centrality(eigenvector_centrality))

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot the top centrality #
#~~~~~~~~~~~~~~~~~~~~~~~~~#
import matplotlib.pyplot as plt
os.chdir("multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/figures/centrality_analysis_TFonly")

def plot_centrality(centrality, title, filename, n=50):
    # Sort the centrality dictionary by values (largest to smallest)
    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)

    # Convert the list to a DataFrame and save to CSV
    df = pd.DataFrame(sorted_centrality, columns=['Node', 'Centrality'])
    df.to_csv(title + ".csv", index=False)
    
    # Limit the number of nodes to the top n
    sorted_centrality = sorted_centrality[:n]

    # Separate the keys (nodes) and values (centralities)
    nodes, centralities = zip(*sorted_centrality)

    # Plot the bar chart
    plt.figure(figsize=(20, 5))
    plt.bar(nodes, centralities, color="skyblue")
    plt.title(title)
    plt.xlabel("Nodes")
    plt.ylabel("Centrality")
    plt.xticks(rotation=90)

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# Call the plot_centrality function for each centrality measure, providing a filename for each plot
plot_centrality(in_degree_centrality, "In-degree Centrality", "in_degree_centrality_top50.png")
plot_centrality(out_degree_centrality, "Out-degree Centrality", "out_degree_centrality_top50.png")
plot_centrality(closeness_centrality, "Closeness Centrality", "closeness_centrality_top50.png")
plot_centrality(betweenness_centrality, "Betweenness Centrality", "betweenness_centrality_top50.png")
plot_centrality(eigenvector_centrality, "Eigenvector Centrality", "eigenvector_centrality_top50.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot all with highlight #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

def plot_centrality_with_highlight(centrality, title, filename, highlighted_nodes=None):
    # Sort the centrality values
    centrality_sorted = {k: v for k, v in sorted(centrality.items(), key=lambda item: item[1], reverse=True)}
    
    # Set the colors for the bars
    colors = ['skyblue']*len(centrality_sorted.keys())
    if highlighted_nodes:
        for node, color in highlighted_nodes.items():
            if node in centrality_sorted.keys():
                colors[list(centrality_sorted.keys()).index(node)] = color
    
    # Create the barplot
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(centrality_sorted)), centrality_sorted.values(), color=colors, width=1)
    plt.xlabel('Nodes')
    plt.ylabel('Centrality')
    plt.title(title)

    # Save the plot to a file
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # Show the plot
    plt.close()
    
# Call the function for each centrality measure and save the plots
highlighted_nodes = {'foxd3': 'orange', 'tfap2a': 'purple'}

plot_centrality_with_highlight(in_degree_centrality, 'In-Degree Centrality', 'in_degree_centrality_with_highlight.png', highlighted_nodes)
plot_centrality_with_highlight(out_degree_centrality, 'Out-Degree Centrality', 'out_degree_centrality_with_highlight.png', highlighted_nodes)
plot_centrality_with_highlight(closeness_centrality, 'Closeness Centrality', 'closeness_centrality_with_highlight.png', highlighted_nodes)
plot_centrality_with_highlight(betweenness_centrality, 'Betweenness Centrality', 'betweenness_centrality_with_highlight.png', highlighted_nodes)
plot_centrality_with_highlight(eigenvector_centrality, 'Eigenvector Centrality', 'eigenvector_centrality_with_highlight.png', highlighted_nodes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot centrality for each clusters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# So we can know the top TFs for each "module"

# Read the df_clusters DataFrame
df_clusters = pd.read_csv("multiome/analysis_newref/GRN_scenicplus/ncall_newmotif/figures/cluster_regulons/clustering_results_gb.csv")
# Get unique values in the 'Agglomerative' column
unique_agglomerative_values = df_clusters['Agglomerative'].unique()
# Iterate over unique values in the 'Agglomerative' column
for value in unique_agglomerative_values:
    # Filter rows based on the current Agglomerative value
    filtered_data = df_clusters[df_clusters['Agglomerative'] == value]

    # Subset the Regulon column
    regulon_subset = filtered_data['Regulon']

    # Split the names by the first '_' and export the first part to a new variable
    tf_list = regulon_subset.str.split('_', n=1, expand=True)[0].tolist()

    # Filter the out-degree centrality by the TFs in the tf_list
    out_degree_centrality_filtered = {k: v for k, v in out_degree_centrality.items() if k in tf_list}

    # Sort the filtered out-degree centrality by values
    out_degree_centrality_sorted = dict(sorted(out_degree_centrality_filtered.items(), 
                                        key=lambda x: x[1], reverse=True))

    # Plot the out-degree centrality for the TFs in the tf_list
    plt.figure(figsize=(10, 5))
    plt.bar(out_degree_centrality_sorted.keys(), out_degree_centrality_sorted.values())
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('TFs')
    plt.ylabel('Out-Degree Centrality')
    plt.title(f'Out-Degree Centrality for TFs in Agglomerative Cluster {value}')
    plt.savefig(f'out_degree_centrality_cluster_{value}_genebased.png', dpi=300, bbox_inches='tight')
    plt.close()

