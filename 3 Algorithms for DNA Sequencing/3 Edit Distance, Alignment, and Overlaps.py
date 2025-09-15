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
