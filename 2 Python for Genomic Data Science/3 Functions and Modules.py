# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

def GC(DNA):
    '''Compute GC percentage of a DNA sequence.'''
    DNA = DNA.upper()
    nbases = DNA.count('N')
    gcpercent = float(DNA.count('C') + DNA.count('G')) * 100/(len(DNA)-nbases)
    return gcpercent

GC('AAAGTNNAGTCC') # Output: 40.0
help(GC)
# print(nbases) doesn't work because it is a local variable, not global



def ContainsStopCodon(DNA, frame=0):
    '''Check if a given DNA sequence contains an in-frame stop codon.'''
    
    stop_codon_found = False
    stop_codons = ['tga', 'tag', 'taa']
    for i in range(frame, len(DNA), 3):
        codon = DNA[i:i+3].lower()
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found

seq = 'atgagcggccggct'
ContainsStopCodon(seq, 0) # False
ContainsStopCodon(frame=1, DNA=seq) # True



def ReverseComplement(DNA):
    '''Return the reverse complement of the DNA string.'''

    def Reverse(DNA):
        return DNA[::-1]
    
    def Complement(DNA):
        base_complement = {'A': 'T', 'C':'G', 'G': 'C', 'T': 'A'}
        letters = list(DNA)
        letters = [base_complement[base] for base in letters]
        return ''.join(letters)
        
    DNA = DNA.upper()
    reverse = Reverse(DNA)
    reverse_complement = Complement(reverse)
    
    return reverse_complement

ReverseComplement(DNA = 'AGTGTGGGGCG') # Output: 'CGCCCCACACT'


# Swap x and y
# -------------
def swap1(x,y):
    x=y
    y=x
    return(x,y)

def swap2(x,y):
    return(y,x)

def swap3(x,y):
    z=x
    x=y
    y=z
    return(x,y)

def swap4(x,y):
    x,y=y,x
    return(x,y)

x, y = 1, 5
swap1(x,y) # 5, 5
swap2(x,y) # 5, 1
swap3(x,y) # 5, 1
swap4(x,y) # 5, 1



def f1(x):
    if (x > 0):
        x = 3*x 
        x = x / 2
    return x

def f2(x):
    if (x > 0):
        x = 3*x
    x = x / 2
    return x

a, b, c, d, e = -1, 0, 1, 3/2, 10
lst = [f1(a), f2(a), f1(b), f2(b), f1(c), f2(c), f1(d), f2(d), f1(e), f2(e)] 
# [-1, -0.5, 0, 0.0, 1.5, 1.5, 2.25, 2.25, 15.0, 15.0]



def function1(length):
    if length > 0:
        print(length)
        function1(length - 1)
def function2(length):
    while length > 0:
        print(length)
        function2(length - 1)

function1(3) # 3 2 1
function2(3) # 3 2 1 1 1 1...


def compute(n, x, y):
    '''Returns x + n*y'''
    if n==0: 
        return x
    # print(n-1, x+y, y)
    return compute(n-1, x+y, y)

compute(2, 3, 4)
# compute(-1, 3, 4) never terminates



def valid_dna1(dna):
    '''Only checks first position.'''
    for c in dna:
        if c in 'acgtACGT':
            return True
        else:
            return False

def valid_dna2(dna):
    '''Always returns True, but c is in acgtACGT.'''
    for c in dna:
        if 'c' in 'acgtACGT':
            return 'True'
        else:
            return 'False'

def valid_dna3(dna):
    '''Only returns result for last position.'''
    
    for c in dna:
        flag = c in 'acgtACGT'
    return flag

def valid_dna4(dna):
    '''Correct function to determine if every position 
    in the string is in 'acgtACGT.'''
    for c in dna:
        if not c in 'acgtACGT':
            return False
    return True

dna1 = 'actgatcgtgcatgCGTAGT'
dna2 = 'natgcttnc'

valid_dna3(dna1) # True
valid_dna3(dna2) # False



L1 = [1,1,2,3,4,5]
L2 = [1,1,3,5,7,9]
L3 = [i for i in set(L1) if i in L2]
print(L3) # [1, 3, 5]
type(L3) # List



def f(mystring):
    print(message)
    print(mystring)
    message="Inside function now!"
    print(message)
message="Outside function!"
f("Test function:") # Error; cannot access global variable message if referenced after



def afunction(a1 = 1, a2):              # Non-default argument follows default argument
def afunction(a1 = 1, a2, a3 = 3):      # Non-default argument follows default argument
def afunction(a1 = 1, a2 = 2, a3 = 3):  # All arguments are given a default
       


import random
def create_dna(n, alphabet='acgt'):
    return ''.join([random.choice(alphabet) for i in range(n)])
dna = create_dna(1000000)
len(dna) # 1000000



def count1(dna, base):
    i = 0
    for c in dna:
        if c == base:
            i += 1 
    return i

def count2(dna, base):
    i = 0 
    for j in range(len(dna)):
        if dna[j] == base:
            i += 1 
    return i 

def count3(dna, base):
    match = [c == base for c in dna]
    return sum(match)

def count4(dna, base):
    return dna.count(base)

def count5(dna, base):
    return len([i for i in range(len(dna)) if dna[i] == base])

def count6(dna,base):
    return sum(c == base for c in dna)

dna, base = 'actgtcgac', 'c'
count1(dna, base) # 3
count2(dna, base) # 3
count3(dna, base) # 3
count4(dna, base) # 3
count5(dna, base) # 3
count6(dna, base) # 3



# ------------------------------------------------------------
# Modules
# ------------------------------------------------------------

export PYTHONPATH=C/user/my_modules:/opt/other_modules
import my_module
# Checks currenty directory --> /user/my_modules --> /opt/other_modules 
        # --> standard library directories --> site-packages


# Check Python version
import sys
print(sys.version)


# Check directory
import sys
sys.path
 # Output: ['C:\\Users\\olivi\\AppData\\Local\\spyder-6\\envs\\spyder-runtime\\python311.zip',
  # 'C:\\Users\\olivi\\AppData\\Local\\spyder-6\\envs\\spyder-runtime\\DLLs', ...]

# Add directory
sys.path.append("C:/Users/olivi/OneDrive/Python/Genomic Data Science/2 Python for Genomic Data Science")
sys.path

# Put related functions into a script to import all at once
import dnautil # in your directory
dir(dnautil) # list all attributes of the module


# Use function from imported module
dna = 'atgagggctaggt'
dnautil.GC(dna) # 53.8


# Import only some functions from the module
from dnautil import GC, ContainsStopCodon


# Determine which function has the shortest run time
import time
dna = 'actgtcgacacgtacgtacgtagctagctacgtacgtacgtacgtacgtacgtacgtacgtagctactgatcgatcgtacgatcgatcgatcgatcgatcgatcgatcgtagctagctacgtacgttgatcgtctacgtactatatacgctgactgatcgatcgctgctgctcgtcgctagctgactgctgctagctgactgatcgatcgactgac'
base = 'c'
functions = [count1, count2, count3, count4, count5, count6]
times = []

for func in functions:
    start = time.perf_counter()
    func(dna, base)
    end = time.perf_counter()
    times.append(end - start)

print("Execution times (in seconds):", times)
fastest_index = times.index(min(times))
print(f"The fastest function is count{fastest_index + 1} with a time of {times[fastest_index]:.6f} seconds.")



# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
# Multiple modules grouped together

# From a package called bioseq:
import bioseq.dnautil
bioseq.dnautil.GC(dna)
# or
from bioseq import dnautil

from bioseq.fasta.fastautil import fastaseqread
