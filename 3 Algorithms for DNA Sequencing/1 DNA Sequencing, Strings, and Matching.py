import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/3 Algorithms for DNA Sequencing/Data")

# ----------------------------------------------
# Basics/Defintions
# ----------------------------------------------


# Length |S|
# ---------------------
s = 'ACGT'
len(s) # |s| = 4

ε = '' # Epsilon is defined as the empty string
len(ε) # |s| = 0


# Offsets
# ---------------------
s = 'ACGT'
s[0] # Leftmost offset = 'A'
s[2] # Offset at position 2 = 'G'


# Concatenation
# ---------------------
s = 'AACC'
t = 'GGTT'
print(s + t)  # 'AACCGGTT'


# Substring
# ---------------------
s = 'AACCGGTT'
s[2:6] # Substring: 'CCGG



# Prefix
# ---------------------
s = 'AACCGGTT'
s[0:6] # Prefix = 'AACCGG'
       # or s[:6]


# Suffix
# ---------------------
s = 'AACCGGTT'
s[4:8] # Suffix = 'GGTT'
s[-4:] # last 4 characters (start at 4 from the end and go to end)
        # or s[4:len(seq)]

# Joining List Elements
# ---------------------
seqs = ['A', 'C', 'G', 'T']
print(''.join(seqs)) # ACGT



# ----------------------------------------------
# Generating Random Strings
# ----------------------------------------------
import random
random.seed(7) # uses random mechanism but fixes the seed (same answers each run)
random.choice('ACGT') # Output: one of the four letters


# Long version
# --------------
seq=''
for _ in range(10): # _ when you're not using the variable (usually i) in the loop
    seq += random.choice('ACGT')
print(seq) # CTAAAGACAA


# As one line
# --------------
seq2 = ''.join([random.choice('ACGT') for _ in range(10)])
print(seq2) # TTACATAACA
seq[1:3] # TA
seq[:3] # TTA
seq[7:] # ACA 
seq[-3] # C


# ----------------------------------------------
# Find Longest Common Prefix
# ----------------------------------------------
def LongestCommonPrefix(s1, s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]

LongestCommonPrefix('ACCATGT', 'ACCAGAC') # Output: 'ACCA


# ----------------------------------------------
# Determine Whether Strings Match
# ----------------------------------------------
def Match(s1, s2):
    if not len(s1) == len(s2):
        return False
    for i in range(len(s1)):
        if not s1[i] == s2[i]:
            return False
        return True
    
Match('ATGCT', 'ATGCT')  # True
Match('ATGCT', 'GTACGT') # False

'ATGCT' == 'GTACGT'   # Built into Python


# ----------------------------------------------
# Reverse Complement
# ----------------------------------------------
def ReverseComplement(s):
    complement = {'A': 'T', 'C':'G', 'G': 'C', 'T': 'A'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

ReverseComplement('AGTGTGGGGCG') # Output: 'CGCCCCACACT'



# ----------------------------------------------
# Downloading and Parsing Genomes
# ----------------------------------------------
# Import FASTA file containing the lambda phage reference genome
import wget
url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa"
wget.download(url, "lambda_virus.fa")

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = ReadGenome('lambda_virus.fa')
genome[:100]
# 'GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACC'
len(genome) # 48502


# Count frequency of each base
# -----------------------------------------------
counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for base in genome:
    counts[base] += 1
print(counts) # {'G': 12820, 'T': 11986, 'A': 12334, 'C': 11362}


import collections
collections.Counter(genome) # Counter({'G': 12820, 'A': 12334, 'T': 11986, 'C': 11362})


# ----------------------------------------------
# FASTQ Files
# ----------------------------------------------

# Format
# ---------------
# line 1: Name
# line 2: Sequence
# line 3: +
# line 4: Base qualities on ASCII scale [Q = -10log10(p)]


# Phred+33: 
# ---------------
# round Q to interger, add 33, convert to character

def QtoPhred33(Q):
    '''Turn Q into Phred+33 ASCII-encoded quality.'''
    return chr(Q+33)

def Phred33toQ(qual):
    '''Turn Phred+33 ASCII-encoded quality into Q = -10log10(p).'''
    return ord(qual)-33 # converts integer to char according to ASCII table


# ----------------------------------------------
# Analyzing Sequencing Read Quality
# ----------------------------------------------
import wget
url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq"
wget.download(url, "SRR835775_1.first1000.fastq")


# Read file
# --------------
def ReadFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

seqs, quals = ReadFastq('SRR835775_1.first1000.fastq')
print(seqs)
print(quals[:5])


def Phred33toQ(qual):
    '''Turn Phred+33 ASCII-encoded quality into Q = -10log10(p).'''
    return ord(qual)-33 # converts integer to char according to ASCII table

Phred33toQ('J') # 41 (less than 1 in 10,000 chance the base is incorrect)


# Create histogram of quality scores
# ------------------------------------
def CreateHist(qualityStrings):
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = Phred33toQ(phred)
            hist[q] += 1
    return hist
h = CreateHist(quals)
print(h)
# [0, 0, 6178, 0, 0, 54, 108, 574, 345, 83, 193, 124, 79, 165, 49, 236, 184, 327, 
# 514, 238, 531, 254, 313, 798, 992, 888, 1396, 1488, 993, 1752, 3387, 4487, 3248, 
# 5476, 8375, 11814, 4243, 7827, 6579, 8179, 9349, 8180, 0, 0, 0, 0, 0, 0, 0, 0]

# %matplotlib inline # for notebook/console
import matplotlib.pyplot as plt
plt.plot(range(len(h)), h)
plt.show()


# ----------------------------------------------
# Analyzing Sequencing Reads by Position
# ----------------------------------------------
def FindGCByPos(reads):
    ''' Find the GC ratio at each position in the read.'''
    
    # Keep track of the number of G/C bases and the total number of bases at each position
    gc = [0] * 100
    totals = [0] * 100
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
            
    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc

gc = FindGCByPos(seqs)
plt.plot(range(len(gc)), gc)
plt.show()


import collections
count = collections.Counter()
for seq in seqs:
    count.update(seq)
    # Counter({'G': 28742, 'C': 28272, 'T': 21836, 'A': 21132, 'N': 18})


# ----------------------------------------------
# Naive Exact Matching
# ----------------------------------------------
import wget
url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa"
wget.download(url, "phix.fa")

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

genome = ReadGenome('phix.fa')

def Naive(p, t):
    '''Returns offset of every match of p within t (location of matches).
    
    For Text T (x = length |T|) and a Pattern P (y = length |P|), there are y-x+1 possible alignments.
    There are at most x(y-x+1) possible character comparisons. Only happens when every char of P matches every character of T.
    There are at least y-x+1 possible comparisons. This happens when first char isn't in the text anywhere.
    In practice, the number is usually closer to the minimum.'''
    
    occurrences = []
    for i in range(len(t) - len(p) + 1): # loop over all alignments; every position where p could start
        match = True
        for j in range(len(p)): # loop over characters
            if t[i+j] != p[j]: # compare characters
                match = False # mismatch; reject alignment
                break
        if match:
          occurrences.append(i) # all chars match; record
    return occurrences

t = 'AGCTTAGATAGC'
p = 'AG'
Naive(p, t) # [0, 5, 9]


# ---------------------------------
# Match artificial reads to genome
# ---------------------------------
import random
def GenerateReads(genome, numReads, readLen):
    ''' Generate reads from random positions in the given genome. '''
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome)-readLen) - 1
        reads.append(genome[start : start+readLen])
    return reads

# Generate 100 reads of length 100
# ----------------------------------
reads = GenerateReads(genome, 100, 100)

# Count how many reads match the genome exactly
# ----------------------------------
numMatched = 0
for r in reads:
    matches = Naive(r, genome)
    if len(matches) > 0:
        numMatched += 1
        
print('%d / %d reads matched the genome exactly!' % (numMatched, len(reads)))
# 100 / 100 reads matched the genome exactly!


# ---------------------------------
# Match real reads to genome
# ---------------------------------
import wget
url = "https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus_2.fa"
wget.download(url, "lambda_virus.fa")

def ReadFastq(filename):
    sequences = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            fh.readline() # skip base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
    return sequences

import collections
phix_reads = ReadFastq('ERR266411_1.first1000.fastq')
count = collections.Counter()
for read in phix_reads:
    count.update(read)
count # Output: Counter({'T': 30531, 'A': 28426, 'C': 21890, 'G': 19147, 'N': 6})


# Number of matching reads
# --------------------------
numMatched = 0
n = 0
for r in phix_reads:
    matches = Naive(r, genome)
    n += 1
    if len(matches) > 0:
        numMatched += 1
        
print('%d / %d reads matched the genome exactly!' % (numMatched, n))
# 7 / 1000 reads matched the genome exactly!
# Maybe sequencing errors?


# Try matching just the first 30 bases of each read
# --------------------------------------------------
numMatched = 0
n = 0
for r in phix_reads:
    r = r[:30]  # just taking the first 30 bases
    matches = Naive(r, genome)
    n += 1
    if len(matches) > 0:
        numMatched += 1
        
print('%d / %d reads matched the genome exactly!' % (numMatched, n))
# 459 / 1000 reads matched the genome exactly!



# Include reverse complement of each read
# ----------------------------------------------
def ReverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

numMatched = 0
n = 0
for r in phix_reads:
    r = r[:30]  # just taking the first 30 bases
    matches = Naive(r, genome)
    matches.extend(Naive(ReverseComplement(r), genome))
    n += 1
    if len(matches) > 0:
        numMatched += 1
        
print('%d / %d reads matched the genome exactly!' % (numMatched, n))
# 932 / 1000 reads matched the genome exactly!


# ----------------------------------------------
# Naive Read Matching with Reverse Complements
# ----------------------------------------------

def NaiveWithRevComp(p, t):
    occurrences = []

    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)

    p_rc = ReverseComplement(p)
    if p_rc != p:  # Avoid double-counting if p is a palindrome
        for i in range(len(t) - len(p_rc) + 1):
            match = True
            for j in range(len(p_rc)):
                if t[i+j] != p_rc[j]:
                    match = False
                    break
            if match:
                occurrences.append(i)

    return occurrences

def ReverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# Examples
# -----------
genome = ReadGenome('lambda_virus.fa')

# How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?
print(len(NaiveWithRevComp("AGGT", genome))) # Output: 306


# How many times does TTAA or its reverse complement occur in the lambda virus genome?
# Hint: TTAA and its reverse complement are equal, so remember not to double count.
print(len(NaiveWithRevComp("TTAA", genome))) # Output: 195


# What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?
print(sorted(NaiveWithRevComp("ACTAAGT", genome))[0]) # Output: 26028


# What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?
print(sorted(NaiveWithRevComp("AGTCGA", genome))[0]) # Output: 450



# ----------------------------------------------
# Naive Approximate Matching
# ----------------------------------------------
def NaiveMismatches(p, t, d):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mismatches += 1
                if mismatches > 2:
                    break  # too many mismatches; reject alignment
        if mismatches <= d:
            occurrences.append(i)
    return occurrences

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# Examples
# ---------
NaiveMismatches(p="ACTTTA", t="ACTTACTTGATAAAGT", d=2) # Output: [0, 4]


p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences = NaiveMismatches(p, t, d=2)
print(occurrences) # Output: [10, 24, 38]


phix_genome = ReadGenome('phix.fa')
print(len(NaiveMismatches('GATTACA', phix_genome, d=2))) # Ouput: 79


# How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches? 
genome = ReadGenome('lambda_virus.fa')
print(len(NaiveMismatches("TTCAAGCC", genome, d=2))) # Output: 191


# What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?
print(sorted(NaiveMismatches("AGGAGGTT", genome, d=2))[0]) # Output: 49



# ----------------------------------------------
# Detect Bad Cycle
# ----------------------------------------------
from collections import defaultdict

def AnalyzeFastqQuality(fastq_path):
    position_qualities = defaultdict(list)
    read_length = None

    with open(fastq_path, 'r') as f:
        line_num = 0
        for line in f:
            line = line.strip()
            line_num += 1

            # Quality scores are every 4th line starting from the 4th
            if line_num % 4 == 0:
                if read_length is None:
                    read_length = len(line)
                for i, char in enumerate(line):
                    position_qualities[i].append(Phred33toQ(char))

    # Compute average quality per cycle
    avg_quality_per_position = {
        pos: sum(qualities) / len(qualities)
        for pos, qualities in position_qualities.items()}

    # Find position with lowest average quality
    worst_cycle = min(avg_quality_per_position, key=avg_quality_per_position.get)

    return worst_cycle

def Phred33toQ(qual):
    '''Turn Phred+33 ASCII-encoded quality into Q = -10log10(p).'''
    return ord(qual)-33 # converts integer to char according to ASCII table

def ReadFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# Example
# ---------
fastq_file = "ERR037900_1.first1000.fastq"
print(AnalyzeFastqQuality(fastq_file)) # Output: 66
