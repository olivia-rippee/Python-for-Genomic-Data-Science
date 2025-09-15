import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/3 Algorithms for DNA Sequencing/Data")


def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# --------------------------------------------------------------------------
# Naive Exact Matching with Number of Character Comparisons and Alignments
# --------------------------------------------------------------------------
def Naive_comparisons(p, t):
    occurrences = []
    comparisons = 0
    alignments = 0

    for i in range(len(t) - len(p) + 1):  # Loop over all alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # Loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # Compare characters
                match = False  # Mismatch; reject alignment
                break
        if match:
            occurrences.append(i)  # All chars match; record
    return occurrences, comparisons, alignments


genome = ReadGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
print(*Naive_comparisons(p, genome)) # [56922] 984143 799954



# --------------------------------------------------------------------------
# Boyer-Moore with Number of Character Comparisons and Alignments
# --------------------------------------------------------------------------
def z_array(s):
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)
    
    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: Zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z

def n_array(s):
    return z_array(s[::-1])[::-1]

def big_l_prime_array(p, n):
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp

def big_l_array(p, lp):
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l

def small_l_prime_array(n):
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

class BoyerMoore(object):
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        
        # Create map from alphabet characters to integers
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
            
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)
    
    def bad_character_rule(self, i, c):
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, i):
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]
    
    def match_skip(self):
        return len(self.small_l_prime) - self.small_l_prime[1]


def boyer_moore_comparisons(p, p_bm, t):
    '''Do Boyer-Moore matching and return occurrences, number of character comparisons, and alignments tried.'''
    i = 0
    occurrences = []
    comparisons = 0
    alignments = 0

    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        alignments += 1  # New alignment attempt

        for j in range(len(p)-1, -1, -1):  # Go backwards
            comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])   # Bad character rule
                skip_gs = p_bm.good_suffix_rule(j)             # Good suffix rule
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift

    return occurrences, comparisons, alignments


# Example
# ---------
genome = ReadGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p_bm = BoyerMoore(p, alphabet='ACGT') # preprocessing
print(*boyer_moore_comparisons(p, p_bm, genome)) # Output: [56922] 165191 127974


# --------------------------------------------------------------------------
# Approximate Matching via Index
# --------------------------------------------------------------------------
from bisect import bisect_left

class Index(object):
    """ Holds a substring index for a text T """
    def __init__(self, t, k):
        self.t = t
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()
    def query(self, p_kmer):
        """ Return index hits for k‐mer = first k characters of p_kmer """
        # Actually, this query treats p_kmer as a full k‐length string
        i = bisect_left(self.index, (p_kmer, -1))
        hits = []
        while i < len(self.index) and self.index[i][0] == p_kmer:
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_match_index(p, t, index, d):
    """
    Find all positions in t where pattern p (length 24) matches with up to d mismatches (substitutions),
    using a k‐mer Index and pigeonhole principle (i.e. partition pattern into n+1 segments).
    
    Returns a tuple (matches, num_index_hits) where
        - matches is a set/list of starting offsets in t
        - num_index_hits is total number of index hits generated when querying the k‐mer index
    """
    assert len(p) == 24
    k = index.k
    assert k == 8  # for 8‐mer index
    seg_count = d + 1
    seg_len = len(p) // seg_count  # since 24 / 3 = 8 when n=2
    matches = set()  # use a set to avoid duplicates
    num_index_hits = 0
    
    for i in range(seg_count):
        start = i * seg_len
        end = start + seg_len
        # In this case, with len = 24, n=2 => seg_len = 8, segments are [0:8], [8:16], [16:24]
        seg = p[start:end]
        seg_matches = index.query(seg)
        num_index_hits += len(seg_matches)
        
        for m in seg_matches:
            # m is the offset in t of an exact match of this segment
            # We imagine aligning full p so that segment [start:end] matches at t[m : m + seg_len]
            # So the candidate full‐pattern alignment's start in t is:
            candidate_start = m - start
            candidate_end = candidate_start + len(p)
            if candidate_start < 0 or candidate_end > len(t):
                continue
            # Count mismatches outside the matching segment
            mismatches = 0
            # before the segment
            for j in range(0, start):
                if p[j] != t[candidate_start + j]:
                    mismatches += 1
                    if mismatches > d:
                        break
            if mismatches > d:
                continue
            # after the segment
            for j in range(end, len(p)):
                if p[j] != t[candidate_start + j]:
                    mismatches += 1
                    if mismatches > d:
                        break
            if mismatches <= d:
                matches.add(candidate_start)
    
    return matches, num_index_hits


# Examples
# ---------
genome = ReadGenome('chr1.GRCh38.excerpt.fasta')
index = Index(genome, 8) # 8-mer index

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
results = approximate_match_index(p, genome, index, d=2)
len(results[0]) # Output: 19

print(results)
# Output: {273669, 681737, 717706, 262042, 635931, 84641, 160162, 724927, 
         # 657496, 160729, 56922, 191452, 551134, 747359, 421221, 147558, 
         # 364263, 465647, 429299} 90


# --------------------------------------------------------------------------
# Approximate Mthcing with Subsequence Index
# --------------------------------------------------------------------------
import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        self.k = k
        self.ival = ival
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i+self.span:ival], i))
        self.index.sort()
    
    def query(self, subseq):
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        
def approximate_search_subseqindex(p, t, subseq_ind, d=0):
    ival = subseq_ind.ival
    span = subseq_ind.span
    occurrences = set()
    num_index_hits = 0

    segment_count = len(p) - span + 1  # slide by 1
    
    for start in range(segment_count):
        subseq = p[start : start + span : ival]
        
        hits = subseq_ind.query(subseq)
        num_index_hits += len(hits)
        
        for hit in hits:
            alignment_start = hit - start
            if alignment_start < 0 or alignment_start + len(p) > len(t):
                continue

            mismatches = 0
            for j in range(len(p)):
                if p[j] != t[alignment_start + j]:
                    mismatches += 1
                    if mismatches > d:
                        break
            
            if mismatches <= d:
                occurrences.add(alignment_start)
    
    return sorted(occurrences), num_index_hits


# Examples
# ---------
ind = SubseqIndex('ATATAT', 3, 2)
print(ind.index) # [('AAA', 0), ('TTT', 1)]

p = 'TTATAT'
print(ind.query(p[0:])) # Output: []
print(ind.query(p[1:])) # Output: [1]


t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = approximate_search_subseqindex(p, t, subseq_ind, d=2)
print(occurrences) # Output: [0, 14]
print(num_index_hits) # Output: 6


# import wget
# url = "http://www.gutenberg.org/ebooks/1110.txt.utf-8"
# wget.download(url, "1110.txt.utf-8")
t = open('1110.txt.utf-8').read()
p = 'English measure backward'
subseq_ind = SubseqIndex(t, 8, 3)

occurrences, num_index_hits = approximate_search_subseqindex(p, t, subseq_ind)
print(occurrences) # [110326]
print(num_index_hits) # 3


p = 'GGCGCGGTGGCTCACGCCTGTAAT'
genome = ReadGenome('chr1.GRCh38.excerpt.fasta')
subseq_ind = SubseqIndex(genome, k=8, ival=3)
occurrences, num_index_hits = approximate_search_subseqindex(p, subseq_ind, index)
print(num_index_hits) # Output: 0
