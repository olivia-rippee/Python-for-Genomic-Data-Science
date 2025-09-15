import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/3 Algorithms for DNA Sequencing/Data")

import wget
url = ""
wget.download(url, "")

# ----------------------------------------------
# Shortest Common Superstring
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

import itertools

def scs(ss):
    '''Returns shortest common superstring of given strings,
       assuming no string is a strict substring of another.'''
        
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
    return shortest_sup

# Example
# -----------
scs(['ACGGATGAGC', 'GAGCGGA', 'GAGCGAG'])
# Output: 'ACGGATGAGCGAGCGGA'


# ----------------------------------------------
# Greedy Shortest Common Superstring
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

import itertools

def scs(ss):
    '''Returns shortest common superstring of given strings,
       assuming no string is a strict substring of another.'''
       
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            #sup += ssperm[i+1][-(len(ssperm[i+1])-olen):]
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def pick_maximal_overlap(reads, k):
    '''Return a pair of reads from the list with a
       maximal suffix/prefix overlap >= k.  Returns
       overlap length 0 if there are no such overlaps.'''
        
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen

def greedy_scs(reads, k):
    '''Greedy shortest-common-superstring merge. Repeat until no edges 
    (overlaps of length >= k) remain.'''
    
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


# Examples
# ----------
greedy_scs(['ABC', 'BCA', 'CAB'], 2) # Output: 'CABCA'

greedy_scs(['ABCD', 'CDBC', 'BCDA'], 1) # Output: 'CDBCABCDA'

scs(['ABCD', 'CDBC', 'BCDA']) # Output: 'ABCDBCDA'


# ----------------------------------------------
# De Bruijn
# ----------------------------------------------

def de_bruijn_ize(st, k):
    '''Return a list holding, for each k-mer, its left k-1-mer and its right k-1-mer in a pair.'''
    edges = []
    nodes = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

nodes, edges = de_bruijn_ize("ACGCGTCG", 3)
nodes 
# Output: {'AC', 'CG', 'GC', 'GT', 'TC'}

edges
# Output: [('AC', 'CG'), ('CG', 'GC'), ('GC', 'CG'), ('CG', 'GT'), ('GT', 'TC'), ('TC', 'CG')]

from graphviz import Digraph

def visualize_de_bruijn(st, k, output_file="de_bruijn_graph"):
    '''Visualize a directed multigraph using graphviz and render it to a file.'''
    nodes, edges = de_bruijn_ize(st, k)

    dot = Digraph(name="DeBruijn_graph", format='png')  # You can change format if needed

    for node in nodes:
        dot.node(node, label=node)
    for src, dst in edges:
        dot.edge(src, dst)

    # Render the graph to a file (output_file.png) and open it
    dot.render(output_file, view=True)
    
# Example
# ----------
visualize_de_bruijn("ACGCGTCG", 3)
