import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/2 Python for Genomic Data Science/Data")

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


# Read sequences
# -------------------------
fasta_sequences = list(SeqIO.parse("dna2.fasta", "fasta"))


# Number of sequences
# -------------------------
len(fasta_sequences)


# Length of sequences
# -------------------------
lengths = {}
for sequence in fasta_sequences:
    lengths[sequence.id] = len(sequence.seq)

ordered = sorted(lengths, key = lengths.get)

for seq_id in ordered:
    print(f"{seq_id}: {lengths[seq_id]} bp")


# Longest ORF
# -------------------------
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

def find_orfs_in_frame(seq, frame_offset):
    '''Find ORFs in a specific forward reading frame (frame_offset = 0,1,2).'''
    seq = seq.upper()
    orfs = []
    i = frame_offset
    while i + 3 <= len(seq):
        codon = seq[i:i+3]
        if codon == START_CODON:
            start = i
            i += 3
            while i + 3 <= len(seq):
                codon = seq[i:i+3]
                if codon in STOP_CODONS:
                    end = i + 3
                    orfs.append({
                        "start": start + 1,            # 1-based
                        "length": end - start,         # in bp
                        "frame": frame_offset + 1,     # 1-based frame
                        "sequence": seq[start:end]
                    })
                    break
                i += 3
        else:
            i += 3
    return orfs

def find_all_forward_orfs(seq_record):
    ''''Find all ORFs in forward frames 1, 2, and 3.'''
    all_orfs = []
    seq = str(seq_record.seq)
    for frame_offset in range(3):  # frames 1, 2, 3
        orfs = find_orfs_in_frame(seq, frame_offset)
        for orf in orfs:
            orf["id"] = seq_record.id
        all_orfs.extend(orfs)
    return all_orfs


# Collect all ORFs
all_orfs = []
orf_by_id = {}

for record in fasta_sequences:
    orfs = find_all_forward_orfs(record)
    orf_by_id[record.id] = orfs
    all_orfs.extend(orfs)


# Longest ORF in reading frame 2
frame2_orfs = [orf for orf in all_orfs if orf["frame"] == 2]
longest_frame2 = max(frame2_orfs, key=lambda x: x["length"], default=None)
if longest_frame2:
    print(f"Q1: Length of longest ORF in reading frame 2: {longest_frame2['length']} bp")
else:
    print("Q1: No ORFs found in reading frame 2.")


# Start position of longest ORF in reading frame 3
frame3_orfs = [orf for orf in all_orfs if orf["frame"] == 3]
longest_frame3 = max(frame3_orfs, key=lambda x: x["length"], default=None)
if longest_frame3:
    print(f"Q2: Starting position of longest ORF in reading frame 3: {longest_frame3['start']}")
else:
    print("Q2: No ORFs found in reading frame 3.")


# Longest ORF in any forward reading frame
longest_overall = max(all_orfs, key=lambda x: x["length"], default=None)
if longest_overall:
    print(f"Q3: Length of longest ORF in any forward reading frame: {longest_overall['length']} bp")
else:
    print("Q3: No ORFs found in any forward reading frame.")


# Longest ORF in sequence 'gi|142022655|gb|EQ086233.1|16'
target_id = "gi|142022655|gb|EQ086233.1|16"
target_orfs = orf_by_id.get(target_id, [])
longest_target = max(target_orfs, key=lambda x: x["length"], default=None)
if longest_target:
    print(f"Q4: Length of longest forward ORF in sequence {target_id}: {longest_target['length']} bp")
else:
    print(f"Q4: No ORFs found in sequence {target_id}.")



# Identify Repeats
# -------------------------
def find_repeats(fasta_sequences, n):
    repeat_counts = defaultdict(int)

    for record in fasta_sequences:
        seq = str(record.seq).upper()
        for i in range(len(seq) - n + 1):
            repeat = seq[i:i + n]
            repeat_counts[repeat] += 1

    return repeat_counts

def most_frequent_repeat(repeat_counts):
    if not repeat_counts:
        return None, 0
    most_common = max(repeat_counts.items(), key=lambda item: item[1])
    return most_common  # (repeat, count)

n = 7
repeats = find_repeats(fasta_sequences, n)
most_common_repeat, count = most_frequent_repeat(repeats)
num_repeats_with_max_count = sum(1 for repeat in repeats.values() if repeat == count)

print(f"Total unique repeats of length {n}: {len(repeats)}")
print(f"Most frequent repeat of length {n}: '{most_common_repeat}' occurs {count} times")
print(f"Number of different repeat of length {n} that occur {count} times: {num_repeats_with_max_count}")

