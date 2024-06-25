## Clustering regulons; region beased
## Author: Zhiyuan, GPT4
## 17 March 2023
## last modified 19 Dec 2023
# reticulate::use_condaenv("cellrank_new")
# !pip install dtaidistance

## regulon_analysis/code/1cluster_regulons_regionBased.py

# This code includes K-means, DBSCAN, Agglomerative Clustering (Hierarchical Clustering), 
# and Hierarchical Clustering with DTW. It saves the cluster figures to the 
# "figures/cluster_regulons" directory and stores the clustering results 
# in a pandas DataFrame, which is then saved as a CSV file named "clustering_results.csv". 
# Note that computing the DTW distance matrix can be time-consuming for large datasets.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from dtaidistance import dtw
from scipy.cluster.hierarchy import dendrogram, linkage
import os
os.chdir("regulon_analysis")

# Create the "figures/cluster_regulons" directory if it doesn't exist
output_dir = "figures/cluster_regulons"
input_dir="GRN_scenicplus/ncall_2023oct_ccb/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
df = pd.read_csv("regulon_analysis/figures/process_regulon_data/eRegulon_auc_rb_filtered.csv")
X=df.to_numpy()
X=X.transpose() 
X.shape
# (186, 16550)

# Standardize the data
k=10

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Dimensionality reduction using PCA
pca = PCA(n_components=50)  # Choose the number of components based on your desired level of dimensionality reduction
X_pca = pca.fit_transform(X_scaled)

# Visualization using t-SNE
tsne = TSNE(n_components=2, random_state = 666)
X_tsne = tsne.fit_transform(X_pca)

df_tsne = pd.DataFrame(X_tsne, columns=['tSNE_1', 'tSNE_2'])
df_tsne['regulon'] =df.columns.tolist()
df_tsne.to_csv(output_dir+"/tSNE_coordinates_rb.csv")


# K-means clustering
kmeans = KMeans(n_clusters=k)  # Choose the number of clusters based on domain knowledge or using the elbow method
kmeans_clusters = kmeans.fit_predict(X_pca)

# DBSCAN clustering
dbscan = DBSCAN(eps=3, min_samples=5)
dbscan_clusters = dbscan.fit_predict(X_pca)

# Agglomerative clustering
agg_clustering = AgglomerativeClustering(n_clusters=k)
agg_clusters_pca = agg_clustering.fit_predict(X_pca)

# Spectral clustering
spectral_clustering = SpectralClustering(n_clusters=k)
spectral_clusters = spectral_clustering.fit_predict(X_pca)

# Agglomerative clustering (Hierarchical Clustering)
agg_clustering = AgglomerativeClustering(n_clusters=k)
agg_clusters_scaled = agg_clustering.fit_predict(X_scaled)

# Dynamic Time Warping (DTW) distance matrix
dtw_distance_matrix = dtw.distance_matrix_fast(X_scaled)
# Hierarchical clustering using DTW distance matrix
hierarchical_dtw = linkage(dtw_distance_matrix, method='complete')
hierarchical_dtw_clusters = AgglomerativeClustering(n_clusters=k, affinity='precomputed', linkage='complete').fit_predict(dtw_distance_matrix)

def plot_clusters_with_legend(title, clusters, index_names, labels=None):
    plt.figure(figsize=(10, 8))
    ax = plt.gca()  # Get the current Axes instance
    unique_clusters = sorted(list(set(clusters)))
    colors = plt.cm.get_cmap('tab10', len(unique_clusters))

    for cluster, color in zip(unique_clusters, colors.colors):
        indices = [i for i, c in enumerate(clusters) if c == cluster]
        label = f"Cluster {cluster}" if labels is None else labels[cluster]
        ax.scatter(X_tsne[indices, 0], X_tsne[indices, 1], c=[color], label=label, s=50)
        
        # Annotate each point with its index name
        for i in indices:
            ax.annotate(index_names[i], 
                        (X_tsne[i, 0], X_tsne[i, 1]), 
                        textcoords="offset points", 
                        xytext=(5,5),  # Offset from each point
                        ha='center',
                        fontsize=6)

    plt.title(title)
    plt.xlabel("t-SNE Component 1")
    plt.ylabel("t-SNE Component 2")
    plt.legend(loc="best")
    plt.savefig(f"{output_dir}/{title}_tSNE_rb.pdf", bbox_inches='tight')
    plt.close()

regulon_names = df.columns.str.split('_', n=2).str[0:2].str.join('_')
plot_clusters_with_legend("K-means Clustering", kmeans_clusters, regulon_names)
plot_clusters_with_legend("DBSCAN Clustering", dbscan_clusters, regulon_names)
plot_clusters_with_legend("Agglomerative Clustering after PCA", agg_clusters_pca, regulon_names)
plot_clusters_with_legend("Agglomerative Clustering", agg_clusters_scaled, regulon_names)
plot_clusters_with_legend("Spectral Clustering", spectral_clusters, regulon_names)
plot_clusters_with_legend("Hierarchical DTW Clustering", hierarchical_dtw_clusters, regulon_names)

# Convert the clustering results into a pandas DataFrame
clustering_results = pd.DataFrame({
    "Regulon": df.columns,
    "K-means": kmeans_clusters,
    "DBSCAN": dbscan_clusters,
    "Agglomerative_pca": agg_clusters_pca,
    "Agglomerative": agg_clusters_scaled,
    "Spectral": spectral_clusters,
    "Hierarchical DTW": hierarchical_dtw_clusters
})

# Save the clustering results to a CSV file
clustering_results.to_csv(f"{output_dir}/clustering_results_rb.csv", index=False)


