# Group data by the comparison file and create a plot for each file
import os
import pandas as pd
import matplotlib.pyplot as plt

# Reload the CSV in case it's modified
file_path = "frequency_similarity.csv"  # Replace with the actual path if different
frequency_data = pd.read_csv(file_path)

# Normalize the Frequency Similarity Scores per Comparison File
frequency_data['Normalized Similarity Score'] = frequency_data.groupby("Comparison File")["Frequency Similarity Score"].transform(
    lambda x: x / x.max()
)

# Create a unified plot for all comparison files
plt.figure(figsize=(12, 8))

# Collect the normalized data for export
normalized_data = []

for comp_file, data in frequency_data.groupby("Comparison File"):
    # Create x-axis labels as a combination of the energy settings
    x_labels = data.apply(
        lambda row: f"D:{row['Diffusion Energy']}, R:{row['Rotation Energy']}, C:{row['Coupling Energy']}", axis=1
    )
    # Plot the normalized similarity scores
    plt.plot(
        x_labels,
        data["Normalized Similarity Score"],
        marker="o",
        linestyle="-",
        label=comp_file,
    )
    # Append normalized data for export
    normalized_data.append(data)

# Customize the plot
plt.xlabel("Simulation Settings (D: Diffusion, R: Rotation, C: Coupling)")
plt.ylabel("Normalized Frequency Similarity Score")
plt.title("Normalized Frequency Similarity Scores for All Comparison Files")
# plt.legend(title="Comparison Files", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

# Save the plot
output_plot = "normalized_frequency_similarity_plot-wo-legend.png"
# plt.savefig(output_plot)
plt.show()

# Combine all normalized data for export
export_df = pd.concat(normalized_data)
output_csv = "normalized_frequency_similarity.csv"
export_df.to_csv(output_csv, index=False)

print(f"Plot saved as '{output_plot}'")
print(f"Normalized data exported as '{output_csv}'")


# Assuming the `frequency_similarity.csv` is available or you can provide the file path.

import seaborn as sns

# Reload the data (if available)
try:
    file_path = "frequency_similarity.csv"  # Replace with the actual path if different
    frequency_data = pd.read_csv(file_path)
    
    # Create a scatter plot for each comparison file
    unique_files = frequency_data["Comparison File"].unique()
    
    for comp_file in unique_files:
        file_data = frequency_data[frequency_data["Comparison File"] == comp_file]
        
        plt.figure(figsize=(8, 6))
        scatter = sns.scatterplot(
            data=file_data,
            x="Diffusion Energy",
            y="Coupling Energy",
            size="Frequency Similarity Score",
            hue="Frequency Similarity Score",
            sizes=(50, 300),
            palette="viridis",
            alpha=0.8
        )
        
        plt.title(f"Frequency Similarity Scatter Plot for {comp_file}")
        plt.xlabel("Diffusion Energy")
        plt.ylabel("Coupling Energy")
        plt.legend(title="Similarity Score", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        
        # Save the plot for each file
        output_file = f"scatter_plot_{comp_file.replace('.csv', '').replace(' ', '_')}.png"
        # plt.savefig(output_file)
        plt.show()
    
    print("Scatter plots saved for each comparison file.")
except FileNotFoundError:
    print("File not found. Please ensure `frequency_similarity.csv` is available.")



