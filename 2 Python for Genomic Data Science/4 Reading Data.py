# ------------------------------------------------------------
# Working with Files
# ------------------------------------------------------------

f=open('myfile', 'r') # read is default

f=open('myfile', 'w') # wwrite is NOT default

f=open('myfile', 'a') # Append to existing file


# File may not exist
# -------------------
try:
    f = open("myfile")
except IOError:
    print("The file does not exist.")


# Read file
# ---------------
for line in f:
    print(line)


f.seek(0) # Go to position 0 in file
f.read() # Doesn't return anything without seek first


f.seek(0)
f.readline() # One line at a time


# Write file
# ---------------
string = 'actagctgctagtgac'
f=open('/Users/olivi/OneDrive/Python/Genomic Data Science/2 Python for Genomic Data Science/test.txt', 'a')
f.write(string)
f.close()


# ------------------------------------------------------------
# FASTA Files
# ------------------------------------------------------------

# Build a dictionary containing all sequences from a FASTA file

try:
    f = open("myfile.fa")
except IOError:
    print("File does not exist.")

seqs = {}
for line in f:
    # discard the newline at the end if any
    line = line.rstrip()
    
    # distinguish header from sequence
    if line[0] =='>': # or line.startswith('>)
        words = line.split()
        name = words[0][1:]
        seqs[name]=''
    else: # sequence, not header
        seqs[name] = seqs[name] + line

f.close


# ------------------------------------------------------------
# Check File Type
# ------------------------------------------------------------
def get_extension1(filename):
    return(filename.split(".")[-1])

def get_extension2(filename):
    import os.path
    return(os.path.splitext(filename)[1])

def get_extension3(filename):
    return filename[filename.rfind('.'):][1:]

filename1 = "myfile.tar.gz"
filename2 = 'myfile'
get_extension1(filename1) # 'gz'     # 'myfile'
get_extension2(filename1) # '.gz'    # ''
get_extension3(filename1) # 'gz'     # ''


# ------------------------------------------------------------
# Retrieving From Dictionaries
# ------------------------------------------------------------
for name, seq in seq.items():
    print(name, seq)

# >python processfasta.py -l 250 myfile.fa # only seqs longer than 250
import sys
print(sys.argv) # ['processfasta.py', 'myfile.fa']


# >python mergefasta.py myfile1.fa myfile2.fa
import sys
tocheck=sys.argv[0] # 'mergefasta.py'
tocheck=sys.argv[1] # 'myfile1.fa'
tocheck=sys.argv[2] # 'myfile2.fa'


# In Home directory with > 'mydir' > files 'foo' and 'bar'
import os
filenames = os.listdir('mydir')
f= open(filenames[0]) # open('foo), but foo is not in current directory


# ------------------------------------------------------------
# Getopt
# ------------------------------------------------------------
#!/user/bin/python
import sys
import getopt

def usage():
    print 
    """processfasta.py: reads a FASTA file and builds a dictionary with all 
    sequences bigger than a given length
    
    processfasta.py [-h] [-l <length>] <filename>
    
    -h              print this message
    -l              filter all seqeunces with a length smaller than <length>
                    (default <length>=0)
    <filename>      the file has to be in FASTA format
    """

o, a = getopt.getopt(sys.argv[1:], 'l:h')
# o = optional; a = required; : means l needs an argument

opts = {}
seqlen = 0;

for k,v in o:
    opts[k] = v
if '-h' in opts.keys():
    usage(); sys.exit
if len(a) < 1:
    usage(); sys.exit("Input FASTA file is missing.")
if '-l' in opts.keys():
    if len(opts['l']) < 0:
        print("Length of sequence should be positive."); sys.exit(0)
    seqlen = opts['l']


# ------------------------------------------------------------
# Standard Streams
# ------------------------------------------------------------
# stdin = standard input
sys.stdin.read() # Ctrl D to end input

# stout = standard output
sys.stdout.write("Some useful output.\n") # Returns string and its length

# stderr = output error messages/diagnostics
sys.stderr.write("Warning: input file was not found\n") # Returns string and its length



# ------------------------------------------------------------
# External Calls
# ------------------------------------------------------------
import subprocess
subprocess.call(["ls", "-l"])
subprocess.call(["tophat", "genome_mouse_idx", "PE_reads_1.fq.gz", "PE_reads_2.fq.gz"])
