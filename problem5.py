import json
from Bio import SeqIO
import csv

# Load Sample A spike protein sequences
samples = list(SeqIO.parse("data/p5-aus-spike-proteins.fasta", "fasta"))

# Load the reference sequence
reference_record = SeqIO.read("data/reference.gb", "genbank")

for feature in reference_record.features:
    if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == "S":
        reference_spike = feature.qualifiers["translation"][0]
        break

# Load variant definitions
with open("data/variants.json", "r") as file:
    variants = json.load(file)

#build dictionary of just nonsynonymous spike mutations
spike_variants = {}
for variant_name, variant_info in variants.items():
    spike_variants[variant_name] = [
        mutation
        for mutation in variant_info["nonsynonymous"]
        if mutation["gene"] == "S"
    ]


# classify each sample into one of 19 variants
results = []
for sample in samples:
    sample_sequence = str(sample.seq)

    # Check whether the sample matches the reference exactly
    if sample_sequence == reference_spike:
        accession = sample.id.split("|")[-1]
        results.append([accession, "NO_MUTATIONS"])
        continue

    # Compare the sample with each variant
    best_variant = "UNKNOWN"
    best_percentage = 0
    for variant_name, mutations in spike_variants.items():

        # Skip variants with no spike mutations
        if len(mutations) == 0:
            continue

        matching_mutations = 0
        known_mutations = 0

        for mutation in mutations:
            position = mutation["pos"] - 1

            # ignore ambiguous amino acids
            if sample_sequence[position] == "X":
                continue

            known_mutations += 1

            if sample_sequence[position] == mutation["right"]:
                matching_mutations += 1

        percentage = matching_mutations / known_mutations * 100

        if percentage > best_percentage:
            best_percentage = percentage
            best_variant = variant_name

    # Require at least 51% matching mutations
    if best_percentage >= 51:
        accession = sample.id.split("|")[-1]
        results.append([accession, best_variant])
    else:
        accession = sample.id.split("|")[-1]
        results.append([accession, "UNKNOWN"])


with open("problem5-classification.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(results)
