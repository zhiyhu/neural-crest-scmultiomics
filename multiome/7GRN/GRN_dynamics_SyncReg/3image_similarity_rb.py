"""
Compute images similarity
Author: Zhiyuan, GPT4
Date: 20 Apr 2023
Last modified: 17 Jan 2024
"""
# reticulate::use_condaenv("keras_env")

## run latent_time_grn/code/process_regulon_data.R before this

# Import necessary libraries
import numpy as np
from tensorflow.keras.preprocessing.image import load_img, img_to_array
from tensorflow.keras.applications.vgg16 import VGG16, preprocess_input
from tensorflow.keras.models import Model
from sklearn.metrics.pairwise import cosine_similarity
import os
import seaborn as sns
import matplotlib.pyplot as plt
os.chdir("latent_time_grn/figures/image_similarity_rb/")

# Function to load images from a given folder
def load_images_from_folder(folder_path):
    images = []
    filenames = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".png"):
            img_path = os.path.join(folder_path, filename)
            img = load_img(img_path, target_size=(224, 224))
            images.append(img)
            filenames.append(os.path.splitext(filename)[0])
    return images, filenames

# Function to preprocess images for VGG16 model
def preprocess_images(images):
    preprocessed_images = []
    for img in images:
        img_array = img_to_array(img)
        img_array = np.expand_dims(img_array, axis=0)
        img_preprocessed = preprocess_input(img_array)
        preprocessed_images.append(img_preprocessed)
    return preprocessed_images

# Function to extract features from preprocessed images using VGG16 model
def extract_features(preprocessed_images, model):
    features = []
    for img in preprocessed_images:
        feature = model.predict(img)
        features.append(feature)
    return features

# Function to calculate similarity between image features
def calculate_similarity(features):
    similarities = cosine_similarity(np.vstack(features))
    return similarities

def add_filenames_to_similarity_matrix(similarities, filenames):
    n = len(filenames)
    similarities_with_filenames = np.empty((n + 1, n + 1), dtype=object)
    
    for i, filename in enumerate(filenames):
        similarities_with_filenames[i + 1, 0] = filename
        similarities_with_filenames[0, i + 1] = filename
    
    similarities_with_filenames[1:, 1:] = similarities
    return similarities_with_filenames


def read_similarity_csv(filepath):
    similarities = np.loadtxt(filepath, delimiter=",")
    return similarities

def plot_similarity_heatmap(similarities, labels):
    # Set font size
    sns.set(font_scale=0.6)

    # Create clustered heatmap
    cluster_map = sns.clustermap(similarities, annot=False, cmap='viridis', cbar=True, square=True,
                                 xticklabels=labels, yticklabels=labels, figsize=(20, 20))

    # Rotate labels
    plt.setp(cluster_map.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(cluster_map.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Set title
    plt.title("Image Similarity Heatmap")

    # Save heatmap to a PDF file
    plt.savefig("clustered_similarity_heatmap.pdf", format='pdf', dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()

# Main function
def main():
    # Path to the image folder
    folder_path = 'latent_time_grn/figures/regulon_auc_rb_scatterplots/'
    images, filenames = load_images_from_folder(folder_path)
    preprocessed_images = preprocess_images(images)

    # Create a VGG16 model and extract features from the images
    base_model = VGG16(weights='imagenet')
    model = Model(inputs=base_model.input, outputs=base_model.get_layer('fc2').output)

    features = extract_features(preprocessed_images, model)
    similarities = calculate_similarity(features)

    # Display similarity scores
    print(similarities)
    
    # Save similarity scores to a CSV file
    np.savetxt("similarities.csv", similarities, delimiter=",")
    
     # Add filenames to the similarity matrix
    similarities_with_filenames = add_filenames_to_similarity_matrix(similarities, filenames)
    
    # Display similarity scores with filenames
    print(similarities_with_filenames)
    
    # Save similarity scores with filenames to a CSV file
    np.savetxt("similarities_with_labels.csv", similarities_with_filenames, delimiter=",", fmt="%s")
    
    # Read in the similarity.csv file
    similarity_filepath = "similarities.csv"
    similarities_from_csv = read_similarity_csv(similarity_filepath)
    print("Similarities from CSV:")
    print(similarities_from_csv)
    # Plot the similarity matrix as a heatmap
    plot_similarity_heatmap(similarities_from_csv, filenames)

# Run the main function
if __name__ == "__main__":
    main()
