import os
os.chdir("C:/Users/olivi/OneDrive/Python/Genomic Data Science/3 Algorithms for DNA Sequencing/Data")

# Boyer-Moore skips alignments that won't be fruitful.
# Try alignments left-to right, and character comparisons right-to-left.
# Bad character rule: Upon mistmatch, skip alignments until either a) mismatch 
    # becomes match or b) P moves past mistmatched character.
# Good suffix rule: Skip until a) there are no mistmaches between P and substring t
    # or b) P moves past t. (t is a substring matched by the inner loop.)

# ----------------------------------------------
# Boyer-Moore
# ----------------------------------------------

import string

def z_array(s):
    '''Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s.'''
    
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
    '''Compile the N array (Gusfield theorem 2.2.2) from the Z array.'''
    return z_array(s[::-1])[::-1]

def big_l_prime_array(p, n):
    '''Compile L' array (Gusfield theorem 2.2.2) using p and N array.
    L'[i] = largest index j less than n such that N[j] = |P[i:]|.'''
    
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp

def big_l_array(p, lp):
    ''' Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| '''
        
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l

def small_l_prime_array(n):
    '''Compile lp' array (Gusfield theorem 2.2.4) using N array.'''
    
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    '''Return tables needed to apply good suffix rule.'''
    
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    ''' Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. '''
    
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    '''Given a full match of P to T, return amount to shift as determined by good suffix rule.'''
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    '''Given pattern string and list with ordered alphabet characters, create and return 
    a dense bad character table. Table is indexed by offsetmthen by character.'''
    
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

class BoyerMoore(object):
    '''Encapsulates pattern and associated Boyer-Moore preprocessing.
    Try alignments left-to right, and character comparisons right-to-left.
    Skips alignmetns that won't be fruitful.
    
    Bad character rule: Upon mistmatch, skip alignments until either a) mismatch 
    becomes match or b) P moves past mistmatched character.
    
    Good suffix rule: Skip until a) there are no mistmaches between P and substring t
    or b) P moves past t. (t is a substring matched by the inner loop.)'''
    
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
        '''Return # skips given by bad character rule at offset i.'''
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, i):
        '''Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule.'''
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]
    
    def match_skip(self):
        '''Return amount to shift in case where P matches T.'''
        return len(self.small_l_prime) - self.small_l_prime[1]

# GCTAGCTCTACGAGTCTA
p = 'TCAA'
p_bm = BoyerMoore(p, alphabet='ACGT') # preprocessing
p_bm.bad_character_rule(2, 'T') # 2 = index of match
# Output: 2

# GCTAGCTCTACGAGTCTA
# ACTA
p = 'ACTA'
p_bm = BoyerMoore(p, alphabet='ACGT') # preprocessing
p_bm.good_suffix_rule(0) # 0 = index of first mismatch
# Output: 3

# ACACGCTCTACGAGTCTA
# ACAC
p = 'ACAC'
p_bm = BoyerMoore(p, alphabet='ACGT') # preprocessing
p_bm.match_skip() # good suffix rule again
# Output: 2

def boyer_moore(p, p_bm, t):
    ''''Do Boyer-Moore matching.'''
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1): # go backwards
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])   # bad character
                skip_gs = p_bm.good_suffix_rule(j)             # good suffix
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


# Examples
# ----------
t = 'GCTAGCTCTACGAGTCTA'
p = 'TCTA'
p_bm = BoyerMoore(p, alphabet='ACGT') # preprocessing

boyer_moore(p, p_bm, t) # Output: [6, 14]


# An algorithm that preprocesses Text is offline, otherwise it is online.
# Naive matching is online (does not preprocessing).
# Boyer-Moore in online (preprocesses Pattern but not Text).
# Web search engines are offline (preprocesses, otherwise it wouldn't be possible in our lifetime).
# An offline is more appropriate for the read alignment problem.
# Indexing is an offline strategy. Given a word/pattern, where can it be found in the text?


# ----------------------------------------------
# Binary Search of the Total k-mer Space
# ----------------------------------------------
# The approximate number of binary-search bisections required to query an 
# ordered index with 32 entries is 5 (since 2^5 = 32).

import bisect
bisect.bisect_left([0, 1, 2, 2, 3, 3, 4], 3) # Where to insert 3 to maintain the sortedness
# Output: 4

# ----------------------------------------------
# Hash Tables
# ----------------------------------------------

t = 'GTGCGTGTGGGG'
table = {'GTG':[0, 4, 6], 'TGC': [1],
         'GCG': [2], 'CGT': [3], 'TGT': [5],
         'TGG': [7], 'GGG': [8, 9, 10]}
table['GGG'] # Output: [8, 9, 10]
table['CGT'] # Output: [3]


# ----------------------------------------------
# Substring Index
# ----------------------------------------------

# Note: subsequence = string of characters also occuring in S in the same order.
# Substrings are strings of consecutive characters.
# Substrings are also subsequences, but subsequences are not necessarily substrings.
# Example: 'AAGT' is a subsequence of 'AACCGGTT' ([0]+[1]+[5]+[7]), but not a substring.

import bisect

class Index(object):
    def __init__(self, t, k):
        '''Create index from all substrings of size 'length'.'''
        
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        '''Return index hits for first k-mer of P.'''
        
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k:i+len(p)]:  # verify that rest of P matches
            offsets.append(i)
    return offsets

# Example
# ---------
t = 'ACTTGGAGATCTTTGAGGCTAGGTATTCGGGATCGAAGCTCATTTCGGGGATCGATTACGATATGGTGGGTATTCGGGA'
p = 'GGTATTCGGGA'
index = Index(t, 4)
print(queryIndex(p, t, index)) # Output: [21, 68]


# ----------------------------------------------
# Approximate Matching Via Boyer-Moore
# ----------------------------------------------
import string

def z_array(s):
    '''Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s.'''
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
    '''Compile the N array (Gusfield theorem 2.2.2) from the Z array.'''
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    '''Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]|'''
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    '''Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]|.'''
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    '''Compile lp' array (Gusfield theorem 2.2.4) using N array.'''
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    '''Return tables needed to apply good suffix rule.'''
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    '''Given a mismatch at offset i, and given L/L' and l' arrays,
       return amount to shift as determined by good suffix rule.'''

    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    '''Given a full match of P to T, return amount to shift as
       determined by good suffix rule.'''
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    '''Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table. Table is indexed by offset
        then by character.'''
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    '''Encapsulates pattern and associated Boyer-Moore preprocessing.'''
    
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
        '''Return # skips given by bad character rule at offset i'''
        
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, i):
        '''Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule.'''
        
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]
    
    def match_skip(self):
        '''Return amount to shift in case where P matches T'''
        return len(self.small_l_prime) - self.small_l_prime[1]

def boyer_moore(p, p_bm, t):
    '''Do Boyer-Moore matching'''
    
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences

def approximate_match(p, t, n):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)

def ReadGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# Examples
# ----------
p = 'AACTTG'
genome = ReadGenome('chr1.GRCh38.excerpt.fasta')
print(approximate_match(p, t, 2)) # Output: [0, 5]
print(approximate_match(p, t, 1)) # Output: [5]
