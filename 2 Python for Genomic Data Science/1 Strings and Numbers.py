# ------------------------------------------------------------
# GC Content
# ------------------------------------------------------------

def GC_Content(DNA):
    '''Calculate the GC percentage of an input string of DNA.
    Input is not case-sensitive.
    
    This wasn't given in the course materials as a function, 
    but I made it one because I couldn't stand not to.'''
    
    DNA = DNA.upper()
    num_c = DNA.count('C')
    num_g = DNA.count('G')
    
    DNA_length = len(DNA)
    
    gc_percent = (num_c + num_g) * 100 / DNA_length
    return(gc_percent)

DNA = 'acgctcgcgcgcgatagctgatcgatcgcgcgcgcttttttttttaaaag'
# or: DNA = raw_input("Enter DNA sequence:")

gc_content = GC_Content(DNA)
print("The DNA sequence's GC content is %5.3f %%" %gc_content) 
                                        # 5 = total digits
                                        # 3 number of digits following decimal
                                        # f: type of value to format
# Output: The DNA sequence's GC content is 52.0 %


# ------------------------------------------------------------
# String Manipulation
# ------------------------------------------------------------
myDNA = 'acgt'
myDNA = myDNA + myDNA
myDNA # Output: 'acgtacgt'
# myDna # Error


DNA="atgctggggact"
DNA[:3]     # Output: 'atg'
DNA         # Output: 'atgctggggact'


DNA='agcagttagcta'
DNA.count('ag') 
# Output: 3


seqlen = '10bp'
seqlen='2'+seqlen
seqlen=seqlen*2
# Output: '210bp210bp'




print('>HSBGPG Human bone gla gene\\transcript "BGP"\nGGCAGATTCCCCCTAGA')

print("""\

>HSBGPG Human bone gla gene\transcript "BGP"

GGCAGATTCCCCCTAGA

""")



dna='agcagttagcta'

# A. Finds the second occurence
o1 = dna.find('atg')
dna.find('atg',o1+1)

# B. Find the last occurence
dna.rfind('atg')

# C. Finds the second occurence
dna.find('atg',dna.find('atg')+1) 




# ------------------------------------------------------------
# Numbers
# ------------------------------------------------------------

17/5 # 8.5 in Python3, 8 in Python2


int(4+6/2+2*2) # Output: 11


lst = [1, 1., 1.0, 1e10, 0x12,"1", "1.0", 100000000000000000, 100000000000000000.0]
for obj in lst:
    print(type(obj))


val1 = 1234567              # 1234567: int
val2 = 1.234567 * 10 ** 6   # 1234567.0: float


a=1
b=2
c=a+b
a = b
a = c
d=a+c
print(a, b, c, d) # Output: 3 2 3 6
