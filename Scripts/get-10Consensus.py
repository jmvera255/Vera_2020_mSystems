from Bio import SeqIO
import numpy as np
import collections
import random
import math
import time
import sys
from defined_functions import *

l = []
seqIDs = []
c = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
	temp = str(record.seq)
	l.append(list(temp))
	n = str(record.id)
	seqIDs.append(n)
	c = c + 1
l = np.array(l)
print(l.shape)
l_align_0 = malign_sharper(l, 3, 83, 98, 1.001, 5, 6)
l_align_0 = malign_iterate_sharper(l_align_0, 1.001, 3, 5, 6)
l_align_0 = malign_iterate_sharp(l_align_0, 1.001, 3, 5, 6)
l_align_0 = malign_iterate(l_align_0, 1.001, 3, 5, 6)
x = 0
print("Alignment# = %s", x)
final_shape = l_align_0.shape[0]
print("# of sequences = ", final_shape)

while True:
	x += 1
	print("Sharp Alignment# = %s", x)
	initial_shape = l_align_0.shape[0]
	(l_align_0, seqIDs) = cherrypick_alignments(l_align_0, l_align_0, 1, 3, 5, 6, seqIDs, sys.argv[2])
	final_shape = l_align_0.shape[0]
	print("# of sequences = ", final_shape)
	l_align_0 = malign_iterate_sharper(l_align_0, 1.001, 3, 5, 6)
	l_align_0 = malign_iterate_sharp(l_align_0, 1.001, 3, 5, 6)
	l_align_0 = malign_iterate(l_align_0, 1.001, 3, 5, 6)
	if initial_shape>final_shape:
		print("Starting New Pass")
	else:
		break

l_align_5 = l_align_0

for x in range(l_align_5.shape[0]):
	for y in range(l_align_5.shape[1]):
		if str(l_align_5[x,y])=='':
			l_align_5[x,y] = '-'
		
while True:
	x += 1
	print(" Alignment# = ", x)
	initial_shape = l_align_0.shape[0]
	(l_align_0, seqIDs) = cherrypick_alignments(l_align_0, l_align_0, 1, 3, 5, 6, seqIDs, sys.argv[2])
	final_shape = l_align_0.shape[0]
	print("# of sequences = ", final_shape)

	print("Info Content = ", seq_info(seq_freq(l_align_0[:, 8:14])))
	if initial_shape>final_shape:
		print("Starting New Pass")
	else:
		break

while True:
	x += 1
	print("Sharper Alignment# = ", x)
	initial_shape = l_align_0.shape[0]
	(l_align_0, seqIDs) = cherrypick_alignments(l_align_0, l_align_0, 1, 3, 5, 6, seqIDs, sys.argv[2])
	final_shape = l_align_0.shape[0]
	print("# of sequences = ", final_shape)

	print("Info Content = ", seq_info(seq_freq(l_align_0[:, 8:14])))
	if initial_shape>final_shape:
		print("Starting New Pass")
	else:
		break

l_align_3 = l_align_0

for x in range(l_align_3.shape[0]):
	for y in range(l_align_3.shape[1]):
		if str(l_align_3[x,y])=='':
			l_align_3[x,y] = '-'
		
out1 = sys.argv[2] + ".-10malign.txt"
with open(out1,"w") as f:
	for x in range(l_align_3.shape[0]):
		for y in range(l_align_3.shape[1]):
			f.write(l_align_3[x,y])
		f.write("\n")

out2 = sys.argv[2] + ".-10malign.seqIDs.txt"
with open(out2,"w") as f2:
	f2.write("\n".join(map(str, seqIDs)))
