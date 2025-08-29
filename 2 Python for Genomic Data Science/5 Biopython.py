import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/2 Python for Genomic Data Science")


# BioPython Module
# ------------------------------------------------------------
import Bio
print(Bio.__version__) # 1.85



# BLAST over Web
# ------------------------------------------------------------
from Bio.Blast import NCBIWWW
fasta_string = open("myseq.fa").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string) # nt = non-redundant nucleotide db
help(NCBIWWW.qblast)



# BLAST Record
# ------------------------------------------------------------
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle) # Nice formatting



# Parsing BLAST Output
# ------------------------------------------------------------
len(blast_record.alignments) # 50 matches (default maximum matches; could be more)

E_VALUE_THRESH = 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
            # print("\n")

# Example Output:
    
    # ****Alignment****
    # sequence: gi|2695885435|gb|OR084852.1| Zaire ebolavirus isolate Ebolavirus/H.sapiens-wt/COD/2020/Mbandaka-RDCEQT000017, partial genome
    # length: 18940
    # e value: 5.21795e-25
    # CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTAAATATATAGTATTGCC
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| || |||||||||||||
    # CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTTAA-ATATAGTATTGCC

    # ****Alignment****
    # sequence: gi|2695885565|gb|OR084865.1| Zaire ebolavirus isolate Ebolavirus/H.sapiens-wt/COD/2020/Mbandaka-RDCEQT000514, partial genome
    # length: 18532
    # e value: 5.21795e-25
    # CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTAAATATATAGTATTGCC
    # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| || |||||||||||||
    # CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTTAA-ATATAGTATTGCC


# Reverse Complement
# ------------------------------------------------------------
from Bio.Seq import Seq

my_seq = Seq("CATGCTACGGTGCTAAAAGCATTA")
print('Reverse complement is %s' % my_seq.reverse_complement())
# Reverse complement is TAATGCTTTTAGCACCGTAGCATG


# Translation
# ------------------------------------------------------------
from Bio.Seq import Seq

DNA = Seq("""TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG""")
RNA = DNA.transcribe()
print(RNA)
# UGGGCCUCAUAUUUAUCCUAUAUACCAUGUUCGUAUGGUGGCGCGAUGUUCUACGUGAAUCCACGUUCGAAGGACAUCAUACCAAAGUCGUACAAUUAGGACCUCGAUAUGGUUUUAUUCUGUUUAUCGUAUCGGAGGUUAUGUUCUUUUUUGCUCUUUUUCGGGCUUCUUCUCAUUCUUCUUUGGCACCUACGGUAGAG


# Translation
# ------------------------------------------------------------
from Bio.Seq import Seq

DNA = Seq("""TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG""")
protein = DNA.translate()
print(protein)
# WASYLSYIPCSYGGAMFYVNPRSKDIIPKSYN*DLDMVLFCLSYRRLCSFLLFFGLLLILLWHLR*



