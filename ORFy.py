# Patrick Ong
# Aleksander Grim

# Project 1: Option 1: Genes and Proteins

# ORFy.py
# This program takes in DNA input from the user
# and generates the following:
### the complementary strand
### the longest gene on each strand
### the mRNA for the longest gene
### the resulting protein

import re

#############

# generates the complementary dna strand, or mRNA, depending on mode
# returns in 5' -> 3' format
def dnaComplement(dna):
	dna = dna.upper()
	dnaComp = ''
	for i in dna:
		if i == 'A':
			dnaComp += 'T'
		elif i == 'T':
			dnaComp += 'A'
		elif i == 'G':
			dnaComp += 'C'
		elif i == 'C':
			dnaComp += 'G'
		else:
			print('\nError: Unexpected character in DNA input !!!\n\n')
			return 'error'

	# reverse string
	dnaRev = dnaComp[::-1]
	return dnaRev


#############

def codonSearch(strands):
	# format:
	# [n][offset]
	# where n is the strand num
	# offset is whether we start on A, T, or G
	# value stored is the length of the gene on that strand starting at specified offset
	#
	# NOTE this assumes no genes tested are 10,000,000 nucleotides long
	geneLen = [[10000000 for x in range(3)] for y in range(2)]
	gene = [['abc' for x in range(3)] for y in range(2)]
	starts = [-1] * 2

	for n in range(2):
		strands[n] = strands[n].upper()

		# NOTE this assumes only one start codon per strand
		start = strands[n].find('ATG')
		starts[n] = start
		if start == -1:
			print('\nError: Start codon not found in one or both DNA strands !!!\n\n')
			return list(-1)

		# note these are empty if the specific stop codon is not found
		stoptaa = [k.start() for k in re.finditer('TAA', strands[n])]
		stoptag = [l.start() for l in re.finditer('TAG', strands[n])]
		stoptga = [m.start() for m in re.finditer('TGA', strands[n])]

		# for all 3 possible starts, A, T, and G
		# find the distance to the first stop codon
		#
		# % 3 to make sure the codon appears as a stop codon
		# rather than something like GCTAAC
		# where the TAA stop is split between neighbors

		for offset in range(3):
			for a in stoptaa:
				if ((a - (start+offset)) % 3 == 0):
					if ((a - start+offset) <= geneLen[n][offset]):
						geneLen[n][offset] = a - (start + offset)
			for b in stoptag:
				if ((b - (start+offset)) % 3 == 0):
					if ((b - (start+offset)) <= geneLen[n][offset]):
						geneLen[n][offset] = b - (start + offset)
			for c in stoptga:
				if ((c - (start+offset)) % 3 == 0):
					if ((c - (start+offset)) <= geneLen[n][offset]):
						geneLen[n][offset] = c - (start + offset)

	# generate the genes themselves

	for x in range(2):
		for y in range(3):
			if ((geneLen[x][y] == 0) | (geneLen[x][y] == 10000000)):
				gene[x][y] = ''
			else:
				gene[x][y] = strands[x][starts[x]+y:starts[x]+geneLen[x][y]+3]

	# find the longest gene on each strand

	ysav = [0]*2
	for x in range(2):
		longest = 0
		for y in range(3):
			if len(gene[x][y]) > longest:
				longest = len(gene[x][y])
				ysav[x] = y
		if x == 0:
			print('The longest gene on the original strand is:\n     5\' ' + gene[x][ysav[x]] + ' 3\'')
		else:
			print('The longest gene on the complementary strand is:\n     5\' ' + gene[x][ysav[x]] + ' 3\'')

	if len(gene[0][ysav[0]]) > len(gene[1][ysav[1]]):
		return gene[0][ysav[0]]
	else:
		return gene[1][ysav[1]]

#############

# returns the single letter amino acid code
# one codon at a time
def codonLookup(codon):
	codon = codon.upper()
	if ((codon == 'UUU') | (codon == 'UUC')):
		return 'F' # phenylalanine
	elif ((codon[0:2] == 'UC') | (codon == 'AGU') | (codon == 'AGC')):
		return 'S' # serine
	elif ((codon == 'UAU') | (codon == 'UAC')):
		return 'Y' # tyrosine
	elif ((codon == 'UAA') | (codon == 'UAG') | (codon == 'UGA')):
		return 'stop' # stop
	elif ((codon == 'UGU') | (codon == 'UGC')):
		return 'C' # cysteine
	elif ((codon == 'UUA') | (codon == 'UUG') | (codon[0:2] == 'CU')):
		return 'L' # leucine
	elif ((codon[0:2] == 'CC')):
		return 'P' # proline
	elif ((codon == 'CAU') | (codon == 'CAC')):
		return 'H' # histidine
	elif ((codon == 'CAA') | (codon == 'CAG')):
		return 'Q' # glutamine
	elif ((codon[0:2] == 'CG') | (codon == 'AGA') | (codon == 'AGG')):
		return 'R' # arginine
	elif ((codon == 'AUU') | (codon == 'AUC') | (codon == 'AUA')):
		return 'I' # isoleucine
	elif ((codon == 'AUG')):
		return 'M' # methionine
	elif ((codon[0:2] == 'AC')):
		return 'T' # threonine
	elif ((codon == 'AAU') | (codon == 'AAC')):
		return 'N' # asparagine
	elif ((codon == 'AAA') | (codon == 'AAG')):
		return 'K' # lysine
	elif ((codon[0:2] == 'GU')):
		return 'V' # valine
	elif ((codon[0:2] == 'GC')):
		return 'A' # alanine
	elif ((codon == 'GAU') | (codon == 'GAC')):
		return 'D' # aspartic acid
	elif ((codon == 'GAA') | (codon == 'GAG')):
		return 'E' # glutamic acid
	elif ((codon[0:2] == 'GG')):
		return 'G' # glycine
	else:
		return 'Z' # UNIDENTIFIED AMINO - ERROR

#############

def RNAtranslate(mRNA):
	amino = codonLookup(mRNA[0:3])
	if (amino == 'stop'):
		return ''
	else:
		return amino + RNAtranslate(mRNA[3:])
	# in a real life setting, we would also check for amino acid Z
	# instead of letting it through, it would throw an error msg

#############

def main():
	# my terminal does not recognize the input string as a string without quotation marks
	# it may be different on other computers
	dna = input('Enter DNA sequence 5\' to 3\', surrounded by quotation marks:\n     ')

	# rudimentary sanitize input
	dna = dna.replace('3','')
	dna = dna.replace('5','')
	dna = dna.replace('\'','')
	dna = dna.replace(' ','') 

	dnaComp = dnaComplement(dna)
	if dnaComp == 'error':
		return

	print("The complementary DNA strand is:\n     5\' " + dnaComp + " 3\'")

	strands = [dna, dnaComp]
	gene = codonSearch(strands)
	if gene[0] == 'error': # error
		return
	
	# generate the mRNA
	mRNA = gene.replace('T','U')

	print('The mRNA for the gene is:\n     5\' ' + mRNA + ' 3\'')
	protein = RNAtranslate(mRNA)

	print('The final protein is:\n     ' + protein)

	return
	
#############

main()


