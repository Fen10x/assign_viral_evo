from Bio import SeqIO
import numpy as np
import pandas as pd

Ebola_records = list(SeqIO.parse("data/p2-ebola-viruses.fasta", "fasta"))
Zika_records = list(SeqIO.parse("data/p2-zika-viruses.fasta", "fasta"))

Ebola_reference = Ebola_records[0]
Zika_reference = Zika_records[0]

#function to create a dataframe from the records
def make_dataframe(samples):
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
    df["collection_date"] = pd.to_datetime(
        df["collection_date"], format="mixed"
    )
    return df

#call the function to create dataframes for Ebola and Zika viruses
Ebola_df = make_dataframe(Ebola_records)
Zika_df = make_dataframe(Zika_records)

#calculate the elapsed days since the first sample for Ebola and Zika viruses
E_start_date = Ebola_df["collection_date"][0]
Z_start_date = Zika_df["collection_date"][0]
Ebola_df["elapsed_days"] = (
    Ebola_df["collection_date"] - E_start_date
).dt.total_seconds() / (24 * 60 * 60)
Zika_df["elapsed_days"] = (
    Zika_df["collection_date"] - Z_start_date
).dt.total_seconds() / (24 * 60 * 60)

#calculate the genetic distances for Ebola and Zika viruses using the Jukes-Cantor correction
from helper_functions import jukes_cantor
E_reference = Ebola_df["sequence"][0]
Z_reference = Zika_df["sequence"][0]

Ebola_df["distance"] = [
    jukes_cantor(E_reference, sequence)
    for sequence in Ebola_df["sequence"]
]

Zika_df["distance"] = [
    jukes_cantor(Z_reference, sequence)
    for sequence in Zika_df["sequence"]
]

#calculate the regression slopes for Ebola and Zika viruses
E_slope, E_intercept = np.polyfit(
    Ebola_df["elapsed_days"][1:],
    Ebola_df["distance"][1:],
    1
)
Z_slope, Z_intercept = np.polyfit(
    Zika_df["elapsed_days"][1:],
    Zika_df["distance"][1:],
    1
)

E_genome_rate = E_slope
Z_genome_rate = Z_slope

print("E genome rate:", E_genome_rate)
print("Z genome rate:", Z_genome_rate)

E_nucleotide_rate = E_slope / len(Ebola_df["sequence"][0])
Z_nucleotide_rate = Z_slope / len(Zika_df["sequence"][0])

print(f"E nucleotide rate: {E_nucleotide_rate:.15f}")
print(f"Z nucleotide rate: {Z_nucleotide_rate:.15f}")


