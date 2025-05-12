"""
import random
import os


def generate_dna_sequence(length, name):
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    if name:
        insert_pos = random.randint(0, length)
        sequence = sequence[:insert_pos] + name + sequence[insert_pos:]
    return sequence, [sequence[i] for i in range(len(sequence)) if sequence[i] in nucleotides]


def calculate_statistics(sequence):
    total = len(sequence)
    if total == 0:
        return 0, 0, 0, 0, 0
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    a_percent = (a_count / total) * 100
    c_percent = (c_count / total) * 100
    g_percent = (g_count / total) * 100
    t_percent = (t_count / total) * 100
    cg_percent = ((c_count + g_count) / total) * 100
    return a_percent, c_percent, g_percent, t_percent, cg_percent


def save_fasta_file(sequence, seq_id, description, filename):
    with open(filename, 'w') as f:
        f.write(f">{seq_id} {description}\n{sequence}\n")


def main():
    # Get user input
    try:
        length = int(input("Enter the sequence length: "))
        if length <= 0:
            print("Sequence length must be positive.")
            return
    except ValueError:
        print("Invalid input. Please enter a valid number.")
        return

    seq_id = input("Enter the sequence ID: ").strip()
    if not seq_id:
        print("Sequence ID cannot be empty.")
        return

    description = input("Provide a description of the sequence: ").strip()
    name = input("Enter your name: ").strip()

    sequence, pure_sequence = generate_dna_sequence(length, name)
    a_percent, c_percent, g_percent, t_percent, cg_percent = calculate_statistics(pure_sequence)

    filename = f"{seq_id}.fasta"
    save_fasta_file(sequence, seq_id, description, filename)

    print(f"\nThe sequence was saved to the file {filename}")
    print("Sequence statistics:")
    print(f"A: {a_percent:.1f}%")
    print(f"C: {c_percent:.1f}%")
    print(f"G: {g_percent:.1f}%")
    print(f"T: {t_percent:.1f}%")
    print(f"%CG: {cg_percent:.1f}%")

    print(f"File will be saved to: {os.path.abspath(filename)}")


if __name__ == "__main__":
    main()
"""

# Program: Random DNA Sequence Generator and Analyzer
# Purpose: Generate a random DNA sequence of user-defined length, insert the user's name,
#          compute nucleotide statistics, save the sequence in FASTA format,
#          plot composition, and detect common restriction sites.
# Context: This script helps teach DNA analysis and quickly checks sequences for GC content, CG/AT ratio, and restriction sites.


import random
import os
import matplotlib.pyplot as plt

# MODIFIED (added dictionary of EcoRI, BamHI, and HindIII recognition motifs to enable automatic restriction site scanning)
RE_SITES = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII":"AAGCTT",
}
# MODIFIED (search the pure DNA sequence for each motif and return their 0-based positions):
def find_restriction_sites(sequence, sites=RE_SITES):
    hits = {}  # initialize result dictionary
    for enzyme, motif in sites.items():  # loop over each enzyme and its motif
        # find all starting positions where motif matches
        positions = [i for i in range(len(sequence) - len(motif) + 1)
                     if sequence[i:i+len(motif)] == motif]
        if positions:  # if any matches found
            hits[enzyme] = positions
    return hits  # return the dict of hits
# JUSTIFICATION:
# Purpose: Let users quickly check if the DNA has common cut sites
# Benefit: Makes cloning easier by showing exactly where these enzymes will cut


def generate_dna_sequence(length, name):
    # choose 'length' nucleotides randomly from A, C, G, T
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    if name:  # if a name was provided
        insert_pos = random.randint(0, length)  # pick a random insertion index
        # insert the name at that position
        sequence = sequence[:insert_pos] + name + sequence[insert_pos:]
    # filter out any non-ACGT characters (i.e. the inserted name)
    pure_sequence = [sequence[i] for i in range(len(sequence)) if sequence[i] in nucleotides]
    return sequence, pure_sequence  # return both versions


def calculate_statistics(sequence):
    total = len(sequence)  # total count of valid nucleotides
    if total == 0:
        return 0, 0, 0, 0, 0, 0  # avoid division by zero
    # count each nucleotide
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    # calculate percentage of each
    a_percent = (a_count / total) * 100
    c_percent = (c_count / total) * 100
    g_percent = (g_count / total) * 100
    t_percent = (t_count / total) * 100
    # GC content as percent
    gc_percent = ((c_count + g_count) / total) * 100
    # ratio of (C+G) to (A+T)
    at_total = a_count + t_count
    cg_at_ratio = (c_count + g_count) / at_total if at_total else 0
    return a_percent, c_percent, g_percent, t_percent, gc_percent, cg_at_ratio


# ORIGINAL:
'''def save_fasta_file(sequence, seq_id, description, filename):
    with open(filename, 'w') as f:
        f.write(f">{seq_id} {description}\n{sequence}\n")'''

# MODIFIED (save sequence to FASTA file with 60 nucleotides per line for better FASTA readability):
def save_fasta_file(sequence, seq_id, description, filename):
    with open(filename, "w") as f:
        # write FASTA header
        f.write(f">{seq_id} {description}\n")
        # write sequence in 60-nt lines
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

# MODIFIED (plot a bar chart of nucleotide percentages for better visualization)
def plot_nucleotide_composition(a, c, g, t, cg_ratio):
    labels = ['A', 'C', 'G', 'T']
    values = [a, c, g, t]
    plt.figure(figsize=(6, 4))  # create figure
    bars = plt.bar(labels, values)  # bar chart
    plt.ylim(0, 100)  # percent scale
    plt.title('Nucleotide Composition (%)')
    plt.ylabel('Percentage')
    # annotate bars with percentage
    for bar, val in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 f"{val:.1f}%", ha='center')
    # annotate CG/AT ratio
    plt.figtext(0.99, 0.01, f"CG/AT ratio: {cg_ratio:.2f}",
                horizontalalignment='right', fontsize=8)
    plt.tight_layout()
    plt.savefig('composition.png')  # save chart file
    plt.close()  # close figure

def main():
    # get desired length from user
    try:
        length = int(input("Enter the sequence length: "))
        if length <= 0:
            print("Sequence length must be positive.")
            return
    except ValueError:
        print("Invalid input. Please enter a valid number.")
        return

    # get FASTA ID
    seq_id = input("Enter the sequence ID: ").strip()
    if not seq_id:
        print("Sequence ID cannot be empty.")
        return

    # get description and user name
    description = input("Provide a description of the sequence: ").strip()
    name = input("Enter your name: ").strip()

    # generate sequence and pure ACGT-only list
    sequence, pure_sequence = generate_dna_sequence(length, name)

    # calculate stats on pure sequence
    a_pct, c_pct, g_pct, t_pct, gc_pct, cg_at = calculate_statistics(pure_sequence)

    # save to FASTA file
    filename = f"{seq_id}.fasta"
    save_fasta_file(sequence, seq_id, description, filename)

    # report file save
    print(f"\nThe sequence was saved to the file {filename}")

    # display statistics
    print("Sequence statistics:")
    print(f"A: {a_pct:.1f}%  C: {c_pct:.1f}%  G: {g_pct:.1f}%  T: {t_pct:.1f}%")
    print(f"GC content: {gc_pct:.1f}%  |  CG/AT ratio: {cg_at:.2f}")

    # find and report restriction sites
    hits = find_restriction_sites(''.join(pure_sequence))
    if hits:
        print("\nRestriction sites found:")
        for enzyme, positions in hits.items():
            # convert to 1-based positions for user-friendliness
            pos_str = ', '.join(str(p+1) for p in positions)
            print(f"  {enzyme} ({RE_SITES[enzyme]}): positions {pos_str}")
    else:
        print("\nNo EcoRI/BamHI/HindIII sites detected in the sequence.")

    # generate and save composition plot
    plot_nucleotide_composition(a_pct, c_pct, g_pct, t_pct, cg_at)
    print("A plot of nucleotide composition was saved as 'composition.png'.")

    # final path report
    print(f"Full path of FASTA file: {os.path.abspath(filename)}")

# run main() if this script executed directly
if __name__ == "__main__":
    main()
