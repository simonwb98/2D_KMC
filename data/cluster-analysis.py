import pandas as pd
import ast
import glob
import numpy as np

def load_and_process_single_csv(file_path):
    """Load and process the single CSV file."""
    data = pd.read_csv(file_path)
    # Convert neighbor frequencies from strings to dictionaries
    data['Averaged Neighbour Frequency'] = data['Averaged Neighbour Frequency'].apply(ast.literal_eval)
    return data

def load_and_process_comparison_csv(file_path):
    """Load and process the comparison CSV file."""
    # Read the metrics and neighbor frequencies
    metrics = pd.read_csv(file_path, nrows=3)
    metrics = metrics.set_index('Metric')['Value'].to_dict()

    frequencies = pd.read_csv(file_path, skiprows=4)
    frequency_dict = dict(zip(frequencies['Degree'], frequencies['Frequency']))

    return metrics, frequency_dict

def calculate_frequency_similarity(sim_freq, comp_freq):
    """Calculate similarity between frequency distributions."""
    degrees = set(sim_freq.keys()).union(comp_freq.keys())
    sim_vector = np.array([sim_freq.get(degree, 0) for degree in degrees])
    comp_vector = np.array([comp_freq.get(degree, 0) for degree in degrees])
    mse = np.mean((sim_vector - comp_vector) ** 2)
    return mse

def calculate_radius_similarity(sim_radius, comp_radius):
    """Calculate similarity for radius of gyration."""
    return abs(sim_radius - comp_radius)

def compare_simulation_to_experiment(single_results, comparison_files):
    """Compare simulation results to experimental data."""
    frequency_scores = []
    radius_scores = []

    for idx, single_row in single_results.iterrows():
        sim_freq = single_row['Averaged Neighbour Frequency']
        sim_radius = single_row['Average Radius of Gyration']

        for file_path in comparison_files:
            comp_metrics, comp_freq = load_and_process_comparison_csv(file_path)
            comp_radius = comp_metrics['Radius of Gyration']

            # Calculate separate similarity scores
            freq_score = calculate_frequency_similarity(sim_freq, comp_freq)
            radius_score = calculate_radius_similarity(sim_radius, comp_radius)

            # Add to frequency similarity table
            frequency_scores.append({
                'Diffusion Energy': single_row['Diffusion Energy'],
                'Rotation Energy': single_row['Rotation Energy'],
                'Coupling Energy': single_row['Coupling Energy'],
                'Comparison File': file_path,
                'Frequency Similarity Score': freq_score
            })

            # Add to radius similarity table
            radius_scores.append({
                'Diffusion Energy': single_row['Diffusion Energy'],
                'Rotation Energy': single_row['Rotation Energy'],
                'Coupling Energy': single_row['Coupling Energy'],
                'Comparison File': file_path,
                'Radius Similarity Score': radius_score
            })

    return pd.DataFrame(frequency_scores), pd.DataFrame(radius_scores)

# File paths
single_file_path = r"C:\Users\User\Desktop\2D_KMC\data\cluster_results.csv"  # Replace with your file path
comparison_files_pattern = r"C:\Users\User\Desktop\research-updates\NetworkQualityAnalysis\*.csv"  # Folder containing comparison files
comparison_file_paths = glob.glob(comparison_files_pattern)

# Load and process files
single_results = load_and_process_single_csv(single_file_path)
frequency_table, radius_table = compare_simulation_to_experiment(single_results, comparison_file_paths)

# Save results to separate files
frequency_table.to_csv("frequency_similarity.csv", index=False)
radius_table.to_csv("radius_similarity.csv", index=False)

print("Frequency similarity results saved to 'frequency_similarity.csv'")
print("Radius similarity results saved to 'radius_similarity.csv'")
