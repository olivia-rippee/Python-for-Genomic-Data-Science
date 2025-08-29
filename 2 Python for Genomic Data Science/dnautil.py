#!usr/bin/python

""""dnautil module contains a few useful functions for dna sequences."""


def GC(DNA):
    '''Compute GC percentage of a DNA sequence.'''
    DNA = DNA.upper()
    nbases = DNA.count('N')
    gcpercent = float(DNA.count('C') + DNA.count('G')) * 100/(len(DNA)-nbases)
    return gcpercent


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