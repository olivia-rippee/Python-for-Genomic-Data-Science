import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/3 Algorithms for DNA Sequencing/Data")

# Hamming distance: minimum substitutions needed to turn one string into the other (same length)

# Edit distance: minimum number of edits (substitutions, insertions, deletions) needed to turn one string into the other
    # edit distance <= hamming distance
    # If x and y are different lengths, edit distance >= length(x) = length(y)
    # edist(αx, βy) = min{ edist(α, β)  + δ(x, y),    [δ(x, y) = 0 if x=y; 0 otherwise]
                         # edist(αx, β) + 1, 
                         # edist(α, βy) + 1 }
                         
# Global alignment: minimize overall edit distance
    # galign(αx, βy) = min{ galign(α, β)  + p(x, y),    [p = penalty matrix]
                          # galign(αx, β) + p(x, -), 
                          # galign(α, βy) + p(-, y) }
                          
# Local alignment: find most similar pair of substrings from x and y
    # lalign(αx, βy) = min{ lalign(α, β)  + s(x, y),    [s = score matrix]
                          # lalign(αx, β) + s(x, -), 
                          # lalign(α, βy) + s(-, y) }


# ----------------------------------------------
# Edit Distance Dynamic Programming
# ----------------------------------------------
def EditDistRecursive(x, y):
    # This implementation is very slow
    if len(x) == 0:
        return len(y)
    elif len(y) == 0:
        return len(x)
    else:
        distHor = EditDistRecursive(x[:-1], y) + 1
        distVer = EditDistRecursive(x, y[:-1]) + 1
        if x[-1] == y[-1]:
            distDiag = EditDistRecursive(x[:-1], y[:-1])
        else:
            distDiag = EditDistRecursive(x[:-1], y[:-1]) + 1
        return min(distHor, distVer, distDiag)

def EditDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]


# Efficiency
# ------------
import time

start_time = time.time()
x = 'shake spea'
y = 'Shakespear'
EditDistRecursive(x, y)
end_time = time.time()
print(f"Execution time: {end_time - start_time:.4f} seconds")
    # Execution time: 4.3450 seconds

start_time = time.time()
x = 'shake spea'
y = 'Shakespear'
EditDistance(x, y)
end_time = time.time()
print(f"Execution time: {end_time - start_time:.4f} seconds")
    # Execution time: 0.0045 seconds


# ----------------------------------------------
# Global Alignment
# ----------------------------------------------
alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8],
         [4, 0, 4, 2, 8],
         [2, 4, 0, 4, 8],
         [4, 2, 4, 0, 8],
         [8, 8, 8, 8, 8]]


# converts from character to its offset in list alphabet
alphabet.index('A') # Output: 0
alphabet.index('G') # Output: 2


# penalty associated with A (from X) mismatching with T (from Y)
score[alphabet.index('A')][alphabet.index('T')] # Output: 4


# penalty associated with C (from X) being deleted in Y
score[alphabet.index('C')][-1] # Output: 8


def globalAlignment(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0] * (len(y)+1))
        
    # Initialize first column
    for i in range(1, len(x)+1):
        D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]

    # Initialize first row
    for j in range(1,len(y)+1):
        D[0][j] = D[0][j-1] + score[-1][alphabet.index(y[j-1])]
        
    # Fill rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
            distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            D[i][j] = min(distHor, distVer, distDiag)
    
    return D[-1][-1]  # return value in bottom right corner


# Example
# ----------
x = 'TATGTCATGC'
y = 'TATGGCAGC'
print(globalAlignment(x,y)) # Output: 12



# First law of assembly: If a suffix of read A is similar to a prefix of read B,
    # then A and B might overlap in the genome.
# Differences can arise from sequencing errors or polyploidy (2 copies).
# Second law of assembly: More coverage leads to more and longer overlaps.


# ----------------------------------------------
# Finding Overlaps
# ----------------------------------------------
def overlap(a, b, min_length=3):
    '''Return length of longest suffix of 'a' matching a prefix of 'b' that is 
    at least 'min_length' characters long. If no such overlap exists, return 0.'''
    
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


# Examples
# ----------
overlap('TTACGT', 'CGTGTGC') # Output: 3

overlap('TTACGT', 'GTGTGC') # Output: 0


# ----------------------------------------------
# Permutations
# ----------------------------------------------
from itertools import permutations
list(permutations([1,2,3], 2))
# [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]


# ----------------------------------------------
# Finding All Overlaps
# ----------------------------------------------
from itertools import permutations

def overlap(a, b, min_length=3):
    '''Return length of longest suffix of 'a' matching a prefix of 'b' that is 
    at least 'min_length' characters long. If no such overlap exists, return 0.'''
    
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def naive_overlap_map(reads, k):
    olaps = {}
    for a, b in permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a, b)] = olen
    return olaps


# Example
# --------
reads = ['ACGGATC', 'GATCAAGT', 'TTCACGGA']
print(naive_overlap_map(reads, 3))
# Output: {('TTCACGGA', 'ACGGATC'): 5, ('ACGGATC', 'GATCAAGT'): 4}



# ----------------------------------------------
# Shortest Edit Distance
# ----------------------------------------------
def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0] * (len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for j in range(len(y)+1):
        D[0][j] = j
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1         # insertion in x / deletion in y
            distVer = D[i-1][j] + 1         # deletion in x / insertion in y
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]      # match
            else:
                distDiag = D[i-1][j-1] + 1  # substitution
            D[i][j] = min(distHor, distVer, distDiag)
    # Bottom-right is exact edit distance (if aligned start to start)
    return D[-1][-1]


def bestMatchEditDistance(p, t):
    """
    Returns the minimal edit distance between pattern p and any substring of t,
    allowing p to slide anywhere in t (fixed-length sliding window).
    """
    m = len(p)
    n = len(t)
    best = None
    for i in range(0, n - m + 1):
        window = t[i : i + m]
        d = editDistance(p, window)
        if (best is None) or (d < best):
            best = d
    return best


# If you want to allow p to match with insertions/deletions in t as well (so that t-substring can be longer/shorter),
# you could do a more general dynamic programming approach similar to approximate matching:
def approxMatchEditDistance(p, t):
    """
    Returns the minimal edit distance allowing insertion/deletion in both p and t,
    so pattern may align partially at any place in t.
    This is analogous to the method shown in “approximate matching” (dynamic programming where we allow
    start anywhere in t).
    """
    m = len(p)
    n = len(t)
    # D[i][j] = edit distance between p[:i] and t[:j]
    D = [[0] * (n + 1) for _ in range(m + 1)]
    # initialize first column: cost of deleting all of p up to i
    for i in range(1, m + 1):
        D[i][0] = i
    # initialize first row: aligning empty p to prefixes of t has 0 cost (so match can start anywhere)
    for j in range(0, n + 1):
        D[0][j] = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if p[i-1] == t[j-1] else 1
            D[i][j] = min(
                D[i-1][j] + 1,        # delete from p
                D[i][j-1] + 1,        # insert into p (or delete from t)
                D[i-1][j-1] + cost    # substitution / match
            )
    # best match is minimum in the last row (i = m), among all j
    best = min( D[m][j] for j in range(0, n + 1) )
    return best

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# Examples
# ---------
p = "GCGTATGC"
t = "TATTGGCTATACGGTT"
print(approxMatchEditDistance(p, t))  # Output: 2


genome = ReadGenome('chr1.GRCh38.excerpt.fasta')

p = 'GCTGATCGATCGTACG'
print(approxMatchEditDistance(p, genome))  # Output: 3

p = 'GATTTACCAGATTGAG'
print(approxMatchEditDistance(p, genome))  # Output: 2


# ----------------------------------------------
# Optimized Overlaps
# ----------------------------------------------
def overlap(a, b, min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_kmer_index(reads, k):
    """Build a k-mer to reads index."""
    kmer_dict = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = set()
            kmer_dict[kmer].add(read)
    return kmer_dict

def overlap_all_pairs(reads, k):
    """Find all (a, b) read pairs with suffix(a) overlapping prefix(b) by ≥ k."""
    kmer_dict = build_kmer_index(reads, k)
    overlaps = set()
    outgoing_reads = set()

    for read in reads:
        suffix = read[-k:]
        candidates = kmer_dict.get(suffix, set())
        for other in candidates:
            if read != other:
                olen = overlap(read, other, min_length=k)
                if olen > 0:
                    overlaps.add((read, other))
                    outgoing_reads.add(read)  # Read has an outgoing edge

    return overlaps, outgoing_reads

def read_fastq(filename):
    sequences = []
    with open(filename) as f:
        while True:
            f.readline()  # Skip name line
            seq = f.readline().strip()  # Sequence line
            f.readline()  # Skip + line
            f.readline()  # Skip quality line
            if not seq:
                break
            sequences.append(seq)
    return sequences


# Examples
# -----------
reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
overlap_all_pairs(reads, 3)
# Output: [('ABCDEFG', 'EFGHIJ'), ('EFGHIJ', 'HIJABC'), ('HIJABC', 'ABCDEFG')]
overlap_all_pairs(reads, 4) #Output:  []


reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
overlap_all_pairs(reads, 4)
# Output: 
    # [('CGTACG', 'TACGTA'), ('CGTACG', 'GTACGT'), ('CGTACG', 'GTACGA'),
    #  ('CGTACG', 'TACGAT'), ('TACGTA', 'ACGTAC'), ('TACGTA', 'CGTACG'),
    #  ('GTACGT', 'TACGTA'), ('GTACGT', 'ACGTAC'), ('ACGTAC', 'GTACGA'),
    #  ('ACGTAC', 'GTACGT'), ('ACGTAC', 'CGTACG'), ('GTACGA', 'TACGAT')]
    #  {'ACGTAC', 'CGTACG', 'GTACGA', 'GTACGT', 'TACGTA'})

overlap_all_pairs(reads, 5)
# Output: 
    # [('CGTACG', 'GTACGT'), ('CGTACG', 'GTACGA'), ('TACGTA', 'ACGTAC'),
    #  ('GTACGT', 'TACGTA'), ('ACGTAC', 'CGTACG'), ('GTACGA', 'TACGAT')]
    #  {'ACGTAC', 'CGTACG', 'GTACGA', 'GTACGT', 'TACGTA'})


reads = read_fastq("ERR266411_1.for_asm.fastq")
overlaps, outgoing_reads = overlap_all_pairs(reads, k=30)
print("Overlapping read pairs (edges):", len(overlaps))           # Output: 904746
print("Reads with outgoing edges (nodes):", len(outgoing_reads))  # Output: 7161
