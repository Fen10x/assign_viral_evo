from helper_functions import jukes_cantor, kimura_two_parameter
from Bio import Entrez, SeqIO
import pandas as pd

"""
Testing the jukes_cantor function with two simple sequences
string1 = "AACT-GGT-GA"
string2 = "AACT-GC-AGA"
print(jukes_cantor(string1, string2))
"""

#extract the collections dates
samples = list(SeqIO.parse("data/p1-sars-cov-2-variants.fasta", "fasta"))
records = []
for sample in samples:
    fields = sample.description.split("|")

    records.append({
        "sample_id": fields[0],
        "country": fields[1],
        "collection_date": fields[2],
        "sequence": str(sample.seq).upper()
    })

df = pd.DataFrame(records)
#convert date to timestamp
df["collection_date"] = pd.to_datetime(
    df["collection_date"],
    format="mixed"
)

# Find the earliest collected sample
sample_A = df.loc[df["collection_date"].idxmin()]
start_date = sample_A["collection_date"]

"""
# code to download th reference genome and save it to a file
Entrez.email = "z5308203@ad.unsw.edu.au"
ref_id = "NC_045512.2"
with Entrez.efetch(
    db="nuccore",
    id=ref_id,
    rettype="gb",
    retmode="text"
) as handle:
    reference = SeqIO.read(handle, "genbank")

SeqIO.write(reference, "data/reference.gb", "genbank")
"""

reference = SeqIO.read("data/reference.gb", "genbank")
# Extract the collection date for the reference genome and conevert it to a timestamp
collection_date = None
for feature in reference.features:
    if feature.type == "source":
        collection_date = feature.qualifiers.get(
            "collection_date", [None]
        )[0]
        if collection_date:
            break
timestamp = pd.to_datetime(collection_date)

import matplotlib.pyplot as plt
import numpy as np

# calculating JC corrected distance and plotting the results
reference_sequence = str(reference.seq).upper()

distances = []
for seq in df["sequence"]:
    distance = jukes_cantor(reference_sequence, seq)
    distances.append(distance)

df["JC_distance"] = distances

#Calculate elapsed time since first sample was collected
df["elapsed_days"] = (
    df["collection_date"] - start_date
).dt.total_seconds() / (24 * 60 * 60)

df["elapsed_years"] = df["elapsed_days"] / 365.25

#sort the dataframe by collection date
df = df.sort_values("collection_date")

# Remove missing or infinite values before regression
plot_df = df[
    np.isfinite(df["elapsed_days"]) &
    np.isfinite(df["JC_distance"])
].copy()

# Fit linear regression using NumPy
slope, intercept = np.polyfit(
    plot_df["elapsed_days"],
    plot_df["JC_distance"],
    1
)

# Calculate R-squared
predicted = slope * plot_df["elapsed_days"] + intercept

ss_res = np.sum((plot_df["JC_distance"] - predicted) ** 2)
ss_tot = np.sum(
    (plot_df["JC_distance"] - plot_df["JC_distance"].mean()) ** 2
)

r_squared = 1 - (ss_res / ss_tot)

# Generate regression line
x_line = np.linspace(
    plot_df["elapsed_days"].min(),
    plot_df["elapsed_days"].max(),
    100
)
y_line = intercept + slope * x_line

# Plot samples and regression line
plt.figure(figsize=(10, 6))

plt.scatter(
    plot_df["elapsed_days"],
    plot_df["JC_distance"],
    alpha=0.7,
    label="Sample genomes"
)

plt.plot(
    x_line,
    y_line,
    color="red",
    label="Linear regression"
)

plt.xlabel("Days since sample A collection")
plt.ylabel("Jukes-Cantor corrected genetic distance")
plt.title("Genetic distance against elapsed time")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("problem1.svg")

print(f"Estimated evolutionary rate: {slope:.6g} substitutions per genome per day")

# Approximate rate per nucleotide per day
genome_length = len(reference_sequence.replace("-", ""))

rate_per_nt_per_day = slope / genome_length

print(
    f"Estimated rate: {rate_per_nt_per_day:.10f} "
    "substitutions per nucleotide per day"
)