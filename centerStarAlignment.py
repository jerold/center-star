#!/usr/bin/python

from sys import stdout
import random

class CellScore:
	match = 0
	mismatch = 1
	insert = 1
	delete = 1

showMatrix = False
showAlignment = False

class Cell:
	def __init__(self, i, j, score):
		self.i = i
		self.j = j
		self.score = score
		self.origin = []
		self.match = False

	def editTypePrimary(self, origin, s1, s2):
		"""Used in the reconstruction of the primary sequence"""
		if origin.i != self.i:
			if origin.j != self.j:
				if self.match:
					return str(s1[origin.i])
				else:
					return str(s1[origin.i])
			else:
				return str(s1[origin.i])
		else:
			return '-'

	def editTypeSecondary(self, origin, s1, s2):
		"""Used in the reconstruction of the secondary sequence"""
		if origin.i != self.i:
			if origin.j != self.j:
				if self.match:
					return str(s2[origin.j])
				else:
					return str(s2[origin.j])
			else:
				return '-'
		else:
			return str(s2[origin.j])

def reconstructAlignment(cell, s1, s2, construct1, construct2):
	"""Recursive method that traces each cell back constructing the alignment based on which cells in the dynamic programming matrix the cell originated"""
	for o in cell.origin:
		if o.i == 0 and o.j == 0:
			if o.i != cell.i and o.j != cell.j:
				if showAlignment:
					print(" c[" + cell.editTypePrimary(o, s1, s2) + construct1 + "]")
					print(" i[" + cell.editTypeSecondary(o, s1, s2) + construct2 + "]")
				return {'i':(cell.editTypePrimary(o, s1, s2) + construct1), 'c':(cell.editTypeSecondary(o, s1, s2) + construct2)}
		else:
			return reconstructAlignment(o, s1, s2, cell.editTypePrimary(o, s1, s2) + construct1, cell.editTypeSecondary(o, s1, s2) + construct2)

def mpa(s, t):
	"""Multiple Pairwise Alignment is used initially to oidentify the center sequence, and then to construct the center star alignment"""
	d = []
	# Init Array
	for i in range(len(s)+1):
		d.append([])
		for j in range(len(t)+1):
			if i == 0:
				d[i].append(Cell(i, j, CellScore.mismatch*j))
			elif j == 0:
				d[i].append(Cell(i, j, CellScore.mismatch*i))
			else:
				d[i].append(Cell(i, j, 0))
	# Score cells and set origins
	for i in range(len(s)):
		for j in range(len(t)):
			# Pick origin cell(s)
			cell = d[i+1][j+1]
			origins = [d[i+1][j], d[i][j], d[i][j+1]]
			originsScores = [d[i+1][j].score + CellScore.insert, d[i][j].score, d[i][j+1].score + CellScore.delete]
			if s[i] == t[j]:
				originsScores[1] = d[i][j].score + CellScore.match
				cell.match = True
			else:
				originsScores[1] = d[i][j].score + CellScore.mismatch
				cell.match = False
			cell.score = min(originsScores)
			for n in range(len(originsScores)):
				if originsScores[n] == cell.score:
					cell.origin.append(origins[n])

	if showMatrix:
		printMatrix(s, t, d)
		print "Score:", d[len(s)][len(t)].score, "\n"

	# From the constructed table we can now produce all optimal alignments
	reconstruction = reconstructAlignment(d[len(s)][len(t)], s, t, '', '')

	# Return both the score of the alignment "D" and the reconstruction
	return {'D':(d[len(s)][len(t)].score), 'reconstruction':reconstruction}

	
def printMatrix(s, t, d):
	"""For easy output viewing we print the alignment matrix"""
	stdout.write("\nMatrix:\n")
	stdout.write("       ")
	for j in range(len(t)):
		stdout.write("   %c" % t[j])
	stdout.write("\n")	
	for i in range(len(s)+1):
		if i > 0:
			stdout.write("   %c" % s[i-1])
		else:
			stdout.write("    ")
		for j in range(len(t)+1):
			stdout.write("%3d " % d[i][j].score)
		stdout.write("  [%d]\n" % i)

if __name__=="__main__":
	# Read in the list of sequences to be aligned
	S = []
	with open('sequences.txt') as f:
		lines = f.readlines()
		for i, line in enumerate(lines):
			# print("[" + str(line.split('\n')[0]) + "]")
			S.append(line.split('\n')[0])

	# Read in the Scoring Matrix, currently the default is PAM scoring
	with open('scoringMatrix.txt') as f:
		lines = f.readlines()
		print("Scoring:")
		for i, line in enumerate(lines):
			option = line.split('\n')[0]
			(cellType, value) = option.split()
			print(str(cellType) + " = " + str(value))
			if cellType == 'match':
				CellScore.match = int(value)
			elif cellType == 'mismatch':
				CellScore.mismatch = int(value)
			elif cellType == 'insert':
				CellScore.insert = int(value)
			elif cellType == 'delete':
				CellScore.delete = int(value)

	# Print all loaded sequences from which we will be finding the Center Star Sequence, and for which we will construct the Center Star Alignment
	print("S = {" + str(", ".join(s for s in S)) + "}")

	showMatrix = True

	# Sequence with the smallest total distance from each other sequence is chosen as our Center Star Sequence
	alignmentMatrix = [[mpa(si, sj) for si in S] for sj in S]
	iSumVector = []
	for i, row in enumerate(alignmentMatrix):
		iSum = 0
		for j, element in enumerate(row):
			if i != j:
				iSum = iSum + element['D']
		iSumVector.append(iSum)

	showMatrix = False

	# Display Alignment Matrix of all alignment scores, and sums
	print("Alignment Matrix:")
	print("[" + "\n[".join(str(", ".join(str(allign['D']) for allign in a)) + "] s" + str(i) + "=" + str(iSumVector[i]) for i, a in enumerate(alignmentMatrix)) + "\n")

	# Produce the Center Star Alignment by adding spaces to previous reconstructions until all sequences are aligned
	for i, iSum in enumerate(iSumVector):
		if iSumVector[i] == min(iSumVector):
			longestSeq = ""
			print("Pairwise Alignments for Sc=s" + str(i) + ", our Center Star Sequence")
			for j, seq in enumerate(S):
				if j != i:
					print("  i=" + str(j))
					print "Sc[" + S[i] + "]"
					print("   " + str(alignmentMatrix[i][j]['reconstruction']['c']))
					print("   " + str(alignmentMatrix[i][j]['reconstruction']['i']))
					if len(alignmentMatrix[i][j]['reconstruction']['i']) > len(longestSeq):
						longestSeq = alignmentMatrix[i][j]['reconstruction']['i']
					print "s" + str(j) + "[" + S[j] + "]\n"

			# Pring the Center Star Alignment
			print("Center Star Alignment:")
			sc = mpa(alignmentMatrix[i][i]['reconstruction']['c'], longestSeq)
			for j, seq in enumerate(S):
				if j == i:
					print("s0: " + str(sc['reconstruction']['c']))
				else:
					print("s" + str(j) + ": " + str(mpa(sc['reconstruction']['c'], alignmentMatrix[i][j]['reconstruction']['i'])['reconstruction']['c']))

