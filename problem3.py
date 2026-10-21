from Bio import SeqIO

records = list(SeqIO.parse("data/p1-sars-cov-2-variants.fasta", "fasta"))
reference = SeqIO.read("data/reference.fasta", "fasta")

alpha_variants  = ["EPI_ISL_3072041","EPI_ISL_1010703","EPI_ISL_1121977","EPI_ISL_1061036","EPI_ISL_1184096",
                   "EPI_ISL_1184504","EPI_ISL_3568432","EPI_ISL_3568427","EPI_ISL_1184097","EPI_ISL_1447289",
                   "EPI_ISL_1615593","EPI_ISL_1660402","EPI_ISL_1660401","EPI_ISL_2404987","EPI_ISL_2462429", 
                   "EPI_ISL_593696"]

delta_variants =  ["EPI_ISL_1904460", "EPI_ISL_1904466", "EPI_ISL_1904469", "EPI_ISL_1904454", "EPI_ISL_1756027", 
                   "EPI_ISL_2404986", "EPI_ISL_2107448", "EPI_ISL_2107529", "EPI_ISL_2650457", "EPI_ISL_2716414", 
                   "EPI_ISL_2633525", "EPI_ISL_2786933", "EPI_ISL_2811809", "EPI_ISL_2695796", "EPI_ISL_2650470", 
                   "EPI_ISL_2811992", "EPI_ISL_2695794", "EPI_ISL_2695799", "EPI_ISL_3281261", "EPI_ISL_3186692",
                   "EPI_ISL_3186471", "EPI_ISL_3186636", "EPI_ISL_3186262", "EPI_ISL_3186487", "EPI_ISL_3280996", 
                   "EPI_ISL_3305513", "EPI_ISL_3072100", "EPI_ISL_6757593", "EPI_ISL_3426420", "EPI_ISL_3539677", 
                   "EPI_ISL_4055434", "EPI_ISL_3834282", "EPI_ISL_3643716", "EPI_ISL_3501929", "EPI_ISL_5535311", 
                   "EPI_ISL_3502476", "EPI_ISL_3502311", "EPI_ISL_3834174", "EPI_ISL_9964573", "EPI_ISL_4747964",
                   "EPI_ISL_4205585", "EPI_ISL_4552276", "EPI_ISL_4205719", "EPI_ISL_4205539", "EPI_ISL_5689721", 
                   "EPI_ISL_4205197", "EPI_ISL_4552130", "EPI_ISL_4747905", "EPI_ISL_5347813", "EPI_ISL_5534640", 
                   "EPI_ISL_5936432", "EPI_ISL_5778590", "EPI_ISL_5534654", "EPI_ISL_5778467", "EPI_ISL_5778492", 
                   "EPI_ISL_6252077", "EPI_ISL_5535134", "EPI_ISL_5305653", "EPI_ISL_7264600", "EPI_ISL_6585779", 
                   "EPI_ISL_6586228", "EPI_ISL_6251927", "EPI_ISL_6585956", "EPI_ISL_7264849", "EPI_ISL_6586173", 
                   "EPI_ISL_6586136", "EPI_ISL_7264591", "EPI_ISL_6757382", "EPI_ISL_7621131", "EPI_ISL_7621166", 
                   "EPI_ISL_8901286", "EPI_ISL_7748984", "EPI_ISL_7621531", "EPI_ISL_7621692", "EPI_ISL_7621124", 
                   "EPI_ISL_7748623", "EPI_ISL_7748687", "EPI_ISL_7621974", "EPI_ISL_8901160", "EPI_ISL_10116719", 
                   "EPI_ISL_8901318", "EPI_ISL_8901337", "EPI_ISL_10116652", "EPI_ISL_8900986", "EPI_ISL_8901124",
                   "EPI_ISL_8901168", "EPI_ISL_8901097", "EPI_ISL_10116931"] 

gene_locations = {
    'S': (21462, 25284),
    'E': (26144, 26372),
    'M': (26422, 27091),
    'N': (28173, 29433)
}

#create dictionary to store the records for the alpha and delta variants
alpha_records = {}
for id in alpha_variants:
    for record in records:
        if id in record.description:
            alpha_records[id] = record
delta_records = {}
for id in delta_variants:
    for record in records:
        if id in record.description:
            delta_records[id] = record

#function to compare sequences and return list of differences
def compare_sequences(seq1, seq2):
    differences = []
    for i, (b1, b2) in enumerate(zip(seq1, seq2)):
        if b1 != b2:
            differences.append((1, b1, b2))
        else:
            differences.append((0, b1, b2))
    return differences

#Calculate the nucleotide mismatches between each sequence and the reference sequence for the variant
alpha_mismatches = []
for id in alpha_variants:
    differences = compare_sequences(str(reference.seq), str(alpha_records[id].seq))
    alpha_mismatches.append((id, differences))

delta_mismatches = []
for id in delta_variants:
    differences = compare_sequences(str(reference.seq), str(delta_records[id].seq))
    delta_mismatches.append((id, differences))

# add up the number of mismatches at each position for the alpha and delta variants
alpha_mismatch_array = [0] * len(alpha_mismatches[0][1])
for i in range(len(alpha_mismatches)):
    for j in range(len(alpha_mismatches[i][1])):
        alpha_mismatch_array[j] += alpha_mismatches[i][1][j][0]

delta_mismatch_array = [0] * len(delta_mismatches[0][1])
for i in range(len(delta_mismatches)):
    for j in range(len(delta_mismatches[i][1])):
        delta_mismatch_array[j] += delta_mismatches[i][1][j][0]

# divide mismatch_array by number of sequences to get the proportion of matches at each position
for i in range(len(alpha_mismatch_array)):
    alpha_mismatch_array[i] /= len(alpha_mismatches)

for i in range(len(delta_mismatch_array)):
    delta_mismatch_array[i] /= len(delta_mismatches)



# find mutations that have occurred in the S gene for either/both variants at greater than 0.5 frequency
alpha_mutations = [0] * len(alpha_mismatch_array)
for i in range(gene_locations['S'][0], gene_locations['S'][1]):
    if alpha_mismatch_array[i] > 0.5:
        alpha_mutations[i] = 1

delta_mutations = [0] * len(delta_mismatch_array)
for i in range(gene_locations['S'][0], gene_locations['S'][1]):
    if delta_mismatch_array[i] > 0.5:
        delta_mutations[i] = 1

"""
print("Alpha variant mutations in S gene at >0.5 frequency:")
for i in range(gene_locations['S'][0], gene_locations['S'][1]):
    if alpha_mutations[i] == 1:
        print(f"mutation at position {i}")
        print(f"mutation of {alpha_mismatches[0][1][i][1]} to {alpha_mismatches[0][1][i][2]}")

print("Delta variant mutations in S gene at >0.5 frequency:")
for i in range(gene_locations['S'][0], gene_locations['S'][1]):
    if delta_mutations[i] == 1:
        print(f"mutation at position {i}")
        print(f"mutation of {delta_mismatches[0][1][i][1]} to {delta_mismatches[0][1][i][2]}")
"""
        
#check if any of the mutations are the same between the two variants
for i in range(gene_locations['S'][0], gene_locations['S'][1]):
    if alpha_mutations[i] == 1 and delta_mutations[i] == 1:
        print(f"mutation at position {i} is the same between the two variants")
        print(f"mutation of {alpha_mismatches[0][1][i][1]} to {alpha_mismatches[0][1][i][2]} in alpha variant")
        print(f"mutation of {delta_mismatches[0][1][i][1]} to {delta_mismatches[0][1][i][2]} in delta variant")






#time to plot!
import matplotlib.pyplot as plt
import numpy as np

fig, (alpha_ax, delta_ax) = plt.subplots(2, 1, figsize=(15, 9))

# X-axis positions
alpha_positions = np.arange(len(alpha_mismatch_array))
delta_positions = np.arange(len(delta_mismatch_array))

# Plot mismatch arrays
alpha_ax.plot(alpha_positions, alpha_mismatch_array)
delta_ax.plot(delta_positions, delta_mismatch_array)

# Gene colours
gene_colors = {
    'S': 'red',
    'E': 'blue',
    'M': 'green',
    'N': 'orange'
}

# Mark gene locations on both subplots
for ax in (alpha_ax, delta_ax):
    for gene, (start, end) in gene_locations.items():
        ax.axvspan(start, end, color=gene_colors[gene], alpha=0.2)
        ax.text(
            (start + end) / 2,
            0.95,
            gene,
            transform=ax.get_xaxis_transform(),
            ha='center',
            va='top'
        )

    ax.set_xlim(20000, len(alpha_mismatch_array))
    ax.set_ylabel('Average mutation frequency')

alpha_ax.set_title('Alpha variant')
delta_ax.set_title('Delta variant')
delta_ax.set_xlabel('Nucleotide position')

plt.tight_layout()
plt.savefig('problem3.svg')
plt.show()
