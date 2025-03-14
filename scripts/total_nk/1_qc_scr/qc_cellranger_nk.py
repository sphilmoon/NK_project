import pandas as pd

# Create a dictionary with the updated data
data = {
    'sample_name': ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10', 's11', 's12'],
    'input_cells': [20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000],
    'estimated_cells': [8583, 16050, 76179, 4291, 26209, 8615, 5924, 9906, 21479, 19460, 18111, 26854],
    'percentage_yield': [42.92, 80.25, 380.90, 21.46, 131.05, 43.08, 29.62, 49.53, 107.40, 97.30, 90.56, 134.27],
    'mean_reads_per_cell': [42724, 28133, 5208, 91431, 18238, 40800, 68532, 34608, 5515, 28051, 18861, 14156],
    'fraction_reads_in_cells': [89.6, 94.3, 98.3, 83.4, 95.4, 76.7, 75.3, 89.0, 81.1, 92.0, 96.1, 89.5],
    'total_genes_detected': [20280, 20679, 23641, 18766, 21598, 22528, 17393, 18394, 20604, 22234, 21206, 22315]
}

# Create a DataFrame
df = pd.DataFrame(data)

# Save to CSV file
df.to_csv('scRNAseq_QC_table_new.csv', index=False)

print("CSV file 'scRNAseq_QC_table_new.csv' has been created with the updated data!")