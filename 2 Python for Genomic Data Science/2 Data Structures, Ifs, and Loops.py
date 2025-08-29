# ------------------------------------------------------------
# Lists []
# ------------------------------------------------------------
# Strings are immutable, Lists are mutable

type([1e-10,(1,2),"BGP",[3]])   
# Ouput: list


grades = [70,80.0,90,100]
(grades[1]+grades[3])/2     
# Output: 90.0


splice_site_pairs = ['GT-AG','GC-AG','AT-AC']
splice_site_pairs[:-1]
# Output: ['GT-AG','GC-AG']

L = [1, 2, 3]
e = 15
L.append(e)
print(L)
# Output: [1, 2, 3, 15]


# ------------------------------------------------------------
# Tuples ()
# ------------------------------------------------------------
# immutable

t = ('a', 'c', 'g', 't')
t.append(('A','C','G','T')) # No append for tuples


# ------------------------------------------------------------
# Sets {}
# ------------------------------------------------------------
# | union; & intersection; - difference

dna='actgtctacgtca'
dna_counts={'t':dna.count('t'),'c':dna.count('c'),'g':dna.count('g'),'a':dna.count('a')}
nt=sorted(dna_counts.keys())
print(nt[-1])
# Output: 't'


L1 = [1,2,3]
L2 = [1,3,5]
L3 = list(set(L1)&set(L2))
print(L3)
# Output: [1, 3] all elements common between L1 and L2 without duplicates


# ------------------------------------------------------------
# Dictionaries {key: value}
# ------------------------------------------------------------

# Accessing an entry
# --------------------
TF_motif = {'SP1': 'gggcgg', 'C/EBP': 'attgccaatt', 'ATF': 'tgacgtca', 
            'c-Myc': 'cacgtg', 'Oct-1': 'atgcaaat'}

print("The recognition sequence for the ATF transcription is %s." % TF_motif['ATF'])
# The recognition sequence for the ATF transcription is tgacgtca.

sorted(TF_motif.keys)     # Access keys
sorted(TF_motif.values)   # Access values


# Adding an entry
# --------------------
TF_motif['AP-1'] = 'tga(g/c)tca'
TF_motif.update({'SP1': 'gggcgg', 'C/EBP': 'attgcgcaat', 'Oct-1': 'atgcaaa'})


# Deleting an entry
# --------------------
del TF_motif['SP1']


dna_counts={'g': 13, 'c': 3, 't': 1, 'a': 16}
del dna_counts['a']
print(dna_counts)
# Output: {'g': 13, 'c': 3, 't': 1}



dna='actgtctacgtca'
dna_counts={'t':dna.count('t'), 'c':dna.count('c'), 'g':dna.count('g'), 'a':dna.count('a')}
max_freq=sorted(dna_counts.values())[-1]
# Output: 4
max_freq_list=sorted(dna_counts.values())
# Output: [2, 3, 4, 4]



someData = { }
someData['cheese'] = 'dairy'
someData['Cheese'] = 'dairy'
someData['Cheese'] = 'Dairy'
someData['cheese'] = 'Dairy'
len(someData) # Output: 2


# ------------------------------------------------------------
# If
# ------------------------------------------------------------

dna='agcgcgggtatatatatgcnccann'

if 'n' or 'N' in dna:
    nbases = dna.count('n') + dna.count('N')
    print('DNA sequence has %d undefined bases.' %nbases)

else:
    print('DNA sequence has no undefined bases.')

# Output: DNA sequence has 3 undefined bases.


motif = 'acgg'
DNA = 'actgctacgacggtgctac'
motif in DNA # True

fold = 2 # condition B
fold = 101 # condition A condition A
fold = 1 # condition A

if fold > 2: 
    print('condition A')
elif fold > 100: 
    print('condition B')
if fold > 2 or fold < 2: 
    print('condition A')
else: 
    print('condition B')


x = 10000000000
if x>10 or x<-10: 
    print('big')
elif x>1000000: 
    print('very big') # Will never print; all x > 1000000 are also > 10
elif x<-1000000: 
    print('very big')
else: 
    print('small')


# ------------------------------------------------------------
# Loops
# ------------------------------------------------------------

# While
# ---------
DNA = 'actgctacgacggtgctac'
pos = DNA.find('gt', 0) # position of donor splice site
while pos > -1:
    print("Donor splice site candidate at position %d" %pos)
    pos=DNA.find('gt', pos+1)
    
    # Donor splice site candidate at position 12

lst = []
i=1
while i<2048:
      i=2*i
      lst.append(i)
print(len(lst)) # Output: 11


seq = 'actg'
i=0
while i<len(seq)+1:
    j=0
    while(j<i+1):
        print(seq[j:i])
        j=j+1
    i=i+1


i = 1
while i < 100:
    if i%2 == 0: 
        break
    i += 1
else:
     i=1000
print(i) # Output: 2


# For
# ---------
motifs={"attccgt", "agggggtttttcg", "gtagc"}
for m in motifs:
    print(m, len(m))
    
    # attccgt 7
    # gtagc 5
    # agggggtttttcg 13


seq = 'actg'
for i in range(len(seq)+1):         # line 1
        for j in range(i):          # line 2
                print(seq[j:i])     # line 3


L1 = [1,2,3,4,5,6,7,8,9]
L2 = [1,3,5,7,9,11,13,15]
L3 = []                              # line 1
for elem in L1:                      # line 2
  if elem in L2 and elem not in L3:  # line 3
   L3.append(elem)                   # line 4   
print(L3)



# mylist=[1,2,2,3,4,5] # True for both
mylist = [1, 2, 2, 3, 4, 5]

d = {}
result = False
for x in mylist:
    if x in d:
        result=True
        break
    d[x] = True
print(result, d) # True {1: True, 2: True}

d = {}
result = False
for x in mylist:
    if not x in d:
        d[x]=True
        continue
    result = True
print(result, d) # True {1: True, 2: True, 3: True, 4: True, 5: True}


# Range
# ---------
for i in range(1, 10, 2): # start, stop, step
    print(i) # 1 3 4 7 9
    

protein = 'SDVIHRYKUUPAKSHGWYVCJRF'
for i in range(len(protein)):
    if protein[i] not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
        print("protein contains invalid amino acid %s at position %d" %(protein[i], i))

    # protein contains invalid amino acid U at position 8
    # protein contains invalid amino acid U at position 9
    # protein contains invalid amino acid J at position 20


protein = 'SDVIHRYKUUPAKSHGWYVCJRF'
for i in range(len(protein)):
    if protein[i] not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
        print("this is not a valid protein sequence.")
        break


protein = 'SDVIHRYKUUPAKSHGWYVCJRF'
corrected_protein=''
for i in range(len(protein)):
    if protein[i] not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
        continue
    corrected_protein=corrected_protein+protein[i]
print("Corrected protein sequence is %s" %corrected_protein)
    
    # Corrected protein sequence is SDVIHRYKPAKSHGWYVCRF


def PrimeNumbers(N):
    for y in range(2, N):
        for x in range(2, y):
            if y % x == 0:
                print(y, 'equals', x, '*', y//x)
                break
        else:
            # loop fell through without finding a factor
            print(y, 'is a prime number')

PrimeNumbers(N=10)
    # 2 is a prime number
    # 3 is a prime number
    # 4 equals 2 * 2
    # 5 is a prime number
    # 6 equals 2 * 3
    # 7 is a prime number
    # 8 equals 2 * 4
    # 9 equals 3 * 3


if motif not in dna:
    pass
else:
    function(motif, dna)


range(1,-23,-3) # 1, -2, -5, -8, -11, -14, -17, -20


# ------------------------------------------------------------
# Boolean Logic and Negation
# ------------------------------------------------------------
# Statement: (rnatype is 'ncRNA' and length>=200) or (rnatype is 'ncRNA' and length==22)
    # --> rnatype is 'ncRNA' and (length>=200 or length==22)
# Negation: rnatype is not 'ncRNA' or (length <200 and length != 22)

