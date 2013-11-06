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

	reconstruction = reconstructAlignment(d[len(s)][len(t)], s, t, '', '')

	return {'D':(d[len(s)][len(t)].score), 'reconstruction':reconstruction}

	
def printMatrix(s, t, d):
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
	S = []
	with open('sequences.txt') as f:
		lines = f.readlines()
		for i, line in enumerate(lines):
			# print("[" + str(line.split('\n')[0]) + "]")
			S.append(line.split('\n')[0])

	print("S = {" + str(", ".join(s for s in S)) + "}")

	showMatrix = True

	alignmentMatrix = [[mpa(si, sj) for si in S] for sj in S]
	iSumVector = []
	for i, row in enumerate(alignmentMatrix):
		iSum = 0
		for j, element in enumerate(row):
			if i != j:
				iSum = iSum + element['D']
		iSumVector.append(iSum)

	showMatrix = False

	print("Alignment Matrix:")
	print("[" + "\n[".join(str(", ".join(str(allign['D']) for allign in a)) + "] s" + str(i) + "=" + str(iSumVector[i]) for i, a in enumerate(alignmentMatrix)) + "\n")


	for i, iSum in enumerate(iSumVector):
		if iSumVector[i] == min(iSumVector):
			longestSeq = ""
			print("Pairwise Alignments for Sc=S[" + str(i) + "]")
			for j, seq in enumerate(S):
				if j != i:
					print("  i=" + str(j))
					print "Sc[" + S[i] + "]"
					print("   " + str(alignmentMatrix[i][j]['reconstruction']['c']))
					print("   " + str(alignmentMatrix[i][j]['reconstruction']['i']))
					if len(alignmentMatrix[i][j]['reconstruction']['i']) > len(longestSeq):
						longestSeq = alignmentMatrix[i][j]['reconstruction']['i']
					# mpa(S[i], S[j])
					print "Si[" + S[j] + "]\n"

			print("Center Star Alignment:")
			sc = mpa(alignmentMatrix[i][i]['reconstruction']['c'], longestSeq)
			for j, seq in enumerate(S):
				if j == i:
					# print("s0: " + str(alignmentMatrix[i][j]['reconstruction']['c']))
					print("s0: " + str(sc['reconstruction']['c']))
				else:
					# print("s" + str(j) + ": " + str(alignmentMatrix[i][j]['reconstruction']['i']))
					print("s" + str(j) + ": " + str(mpa(sc['reconstruction']['c'], alignmentMatrix[i][j]['reconstruction']['i'])['reconstruction']['c']))

