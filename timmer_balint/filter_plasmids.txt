import pandas as pd

# Define file paths
input_file = r"D:\Egyetem\bioinfo\szovegbanyaszat\klebsiella_plasmids_table.tsv"
output_file = r"D:\Egyetem\bioinfo\szovegbanyaszat\filtered_klebsiella_plasmids.tsv"

# Define the resistance and virulence gene columns
resistance_cols = ["AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired", "Gly_acquired", "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired", "Tgc_acquired", "Tmt_acquired", "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired", "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired"]
virulence_cols = ["iucA", "iucB", "iucC", "iucD", "iutA", "iroB", "iroC", "iroD", "iroN", "rmpA", "rmpD", "rmpC", "rmpA2"]

# Load the TSV file
df = pd.read_csv(input_file, sep='\t')

# Filter rows where at least one resistance gene is present AND at least one virulence gene is present
filtered_df = df[(df[resistance_cols] != "-").any(axis=1) & (df[virulence_cols] != "-").any(axis=1)]

# Save the filtered data to a new TSV file
filtered_df.to_csv(output_file, sep='\t', index=False)

print(f"Filtered data saved to {output_file}")
