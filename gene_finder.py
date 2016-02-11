# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Mikhaela

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
	"""Shuffles the characters in the input string
		NOTE: this is a helper function, you do not
		have to modify this in any way """
	return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide): 
	""" Returns the complementary nucleotide
		nucleotide: a nucleotide (A, C, G, or T) represented as a string
		returns: the complementary nucleotide
	>>> get_complement('A')
	'T'
	>>> get_complement('C')
	'G'
	>>> get_complement('T')
	'A'
	get_complement('Q') #this should result in an error 
	
	"""	
	if nucleotide == 'A':
		return 'T' 
	elif nucleotide == 'C':
		return 'G'
	elif nucleotide == 'T': 
		return 'A'
	elif nucleotide == 'G':
		return 'C'
	else: 
		raise ValueError('bad letter')
	
def get_reverse_complement(dna):
	""" Computes the reverse complementary sequence of DNA for the specfied DNA
		sequence

		dna: a DNA sequence represented as a string
		returns: the reverse complementary DNA sequence represented as a string
	>>> get_reverse_complement("ATGCCCGCTTT")
	'AAAGCGGGCAT'
	>>> get_reverse_complement("CCGCGTTCA")
	'TGAACGCGG'
	"""     
	reverse_dna = dna[::-1] #[::-1] reversing a list
	complement = ""
	for nucleotide in reverse_dna: #evaluates along length of dna strand with cureent nucleotide 
		complement += get_complement(nucleotide) #constantly appending (+=) to get complement of nucleotide
	return complement 


def rest_of_ORF(dna):
	""" Takes a DNA sequence that is assumed to begin with a start
		codon and returns the sequence up to but not including the
		first in frame stop codon.  If there is no in frame stop codon,
		returns the whole string. stop codons: TAG TGA TAA 

		dna: a DNA sequence
		returns: the open reading frame represented as a string
	>>> rest_of_ORF("ATGTGAA")
	'ATG'
	>>> rest_of_ORF("ATGAGATAGG")
	'ATGAGA'
	>>> rest_of_ORF("ATGCATGAATGA")
	'ATGCATGAA'
	"""
	stop_codons = ["TAG", "TGA", "TAA"] #what are the stop codons 
	for index in range(0,len(dna),3): #pointing at the begining of each set of 3 in dna string
		codon = dna[index:index+3] #evaluate sets of 3 neuclotides
		if codon in stop_codons:  #if dna has a stop codon
			return dna[0:index] #return the dna string up to the stop codon
	return dna #otherwise, just return the whole dna string


def find_all_ORFs_oneframe(dna):
	""" Finds all non-nested open reading frames in the given DNA
		sequence and returns them as a list.  This function should
		only find ORFs that are in the default frame of the sequence
		(i.e. they start on indices that are multiples of 3).
		By non-nested we mean that if an ORF occurs entirely within
		another ORF, it should not be included in the returned list of ORFs.

		dna: a DNA sequence
		returns: a list of non-nested ORFs
	>>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
	['ATGCATGAATGTAGA', 'ATGTGCCC'] 

	While I was coding I manually tired some other sequences to see if they would work, as I was building the fucntion 

	"""
	start_codon = "ATG"
	index = 0
	result = [] 
	while index < len(dna): #so long as we are within the dna string
		condon = dna[index:index+3] #evaluate codons in sets of three starting at the very begining
		if condon == start_codon : #and if the codon we are looking at is a start codon
			result.append(rest_of_ORF(dna[index:])) #then append to a list of strings from whereever that start codon was until rest_of_ORF finds an end codon
			index = index + len(rest_of_ORF(dna[index:])) #after end codon index from where we left off
		else: #otherwise, if no start codon
			index = index + 3 #add another 3 and go to next codon to look for the start 
	return result
	

def find_all_ORFs(dna):
	""" Finds all non-nested open reading frames in the given DNA sequence in
		all 3 possible frames and returns them as a list.  By non-nested we
		mean that if an ORF occurs entirely within another ORF and they are
		both in the same frame, it should not be included in the returned list
		of ORFs.

		dna: a DNA sequence
		returns: a list of non-nested ORFs

	>>> find_all_ORFs("ATGCATGAATGTAG")
	['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

	if the previous functions worked, and were tested properly, adding extra doctest would not be beneficial 
	"""
	result =  find_all_ORFs_oneframe(dna[0:]) #find_all_ORF_oneframe returns a list
	result += find_all_ORFs_oneframe(dna[1:]) #this is adding a list to a list 
	result += find_all_ORFs_oneframe(dna[2:]) #so that I can see the three possible reading frams 
	return result 

def find_all_ORFs_both_strands(dna):
	""" Finds all non-nested open reading frames in the given DNA sequence on both
		strands.

		dna: a DNA sequence
		returns: a list of non-nested ORFs
	>>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
	['ATGCGAATG', 'ATGCTACATTCGCAT']

	same as find_all_ORFs, the functions this function relys on already have sufficient doctests and I don't think this one needs more tests 
	"""
	result =  find_all_ORFs(dna)
	result += find_all_ORFs(get_reverse_complement(dna))
	return result #now I have a list of fowards and backwards in all reading frames


def longest_ORF(dna):
	""" Finds the longest ORF on both strands of the specified DNA and returns it
		as a string
	>>> longest_ORF("ATGCGAATGTAGCATCAAA")
	'ATGCTACATTCGCAT'
	longest_ORF relies on fucntions that have already been tested 
	"""
	most_long_ORF = ""
	all_ORFs = find_all_ORFs_both_strands(dna) #a list of all ORFs
	for ORF in all_ORFs :
		if len(ORF) > len(most_long_ORF):
			most_long_ORF = ORF 
	return most_long_ORF



def longest_ORF_noncoding(dna, num_trials):
	""" Computes the maximum length of the longest ORF over num_trials shuffles
		of the specfied DNA sequence

		dna: a DNA sequence
		num_trials: the number of random shuffles
		returns: the maximum length longest ORF """
	index = 0
	longest_length = 0
	while index < num_trials:
		length_dna = len(longest_ORF(shuffle_string(dna)))
		if longest_length < length_dna:
			longest_length = length_dna
		index = index + 1 
	return longest_length


def coding_strand_to_AA(dna):
	""" Computes the Protein encoded by a sequence of DNA.  This function
		does not check for start and stop codons (it assumes that the input
		DNA sequence represents an protein coding region).

		dna: a DNA sequence represented as a string
		returns: a string containing the sequence of amino acids encoded by the
				 the input DNA fragment

		>>> coding_strand_to_AA("ATGCGA")
		'MR'
		>>> coding_strand_to_AA("ATGCCCGCTTT")
		'MPA'
		>>> coding_strand_to_AA("TTCTTTTTA")
		'FFL'
	"""
	amino_acid_chain = ""
	index = 0
	while index < len(dna)-2: 
		codon = dna[index:index+3]
		amino_acid_chain += aa_table[codon]
		index = index + 3
	return amino_acid_chain

print coding_strand_to_AA("ATGCCCGCTTT")

def gene_finder(dna):
	""" Returns the amino acid sequences that are likely coded by the specified dna

		dna: a DNA sequence
		returns: a list of all amino acid sequences coded by the sequence dna.
	"""
	threshold = longest_ORF_noncoding(dna,1500)
	all_ORFs = find_all_ORFs_both_strands(dna) # list 
	IF_amino_acid = []
	for ORF in all_ORFs:
		if len(ORF) > threshold: 
			IF_amino_acid.append(coding_strand_to_AA(ORF))
	return IF_amino_acid


from load import load_seq
dna = load_seq("./data/X73525.fa")

#print gene_finder(dna)

if __name__ == "__main__":
	import doctest
	doctest.testmod()
	doctest.run_docstring_examples(coding_strand_to_AA, globals(),verbose = True)

