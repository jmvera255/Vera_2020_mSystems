from __future__ import division
from Bio import SeqIO
import numpy as np
import collections
import random
import math
import time
import csv

def load_fasta_promoters(filename):
	l = []
	records = list(SeqIO.parse(filename, "fasta"))
	for record in records:
		temp = str(record.seq)
		l.append(list(temp))
	l = np.array(l)
	return l

def seq_freq_orig(seq_freq_input):
	freq = []
	for colnbr in range(np.shape(seq_freq_input)[1]):
		col = seq_freq_input[:,colnbr]
		x = collections.Counter(col)
		tot = x['A'] + x['T'] + x['G'] + x['C']
		if tot == 0: tot = 1
		freq.append([float(x['A'])/tot, float(x['T'])/tot, float(x['G'])/tot, float(x['C'])/tot])
	#print('\n'.join(' '.join(map(str,sfreq)) for sfreq in freq))
	return np.transpose(np.array(freq))

def seq_freq(tempmat):
#This one adds random DNA to gaps to avoid any biases during sequence information calculations
	freq = []
	seq_freq_input = tempmat
	for x in range(np.shape(seq_freq_input)[0]):
		for y in range(np.shape(seq_freq_input)[1]):
			if seq_freq_input[x, y] == '':
				n = random.randint(46, 49)
				if n==46:
					seq_freq_input[x, y] = 'A'
				elif n==47:
					seq_freq_input[x, y] = 'T'
				elif n==48:
					seq_freq_input[x, y] = 'G'
				else:
					seq_freq_input[x, y] = 'C'
	for colnbr in range(np.shape(seq_freq_input)[1]):
		col = seq_freq_input[:,colnbr]
		x = collections.Counter(col)
		tot = x['A'] + x['T'] + x['G'] + x['C']
		if tot == 0: tot = 1
		freq.append([float(x['A'])/tot, float(x['T'])/tot, float(x['G'])/tot, float(x['C'])/tot])
	return np.transpose(np.array(freq))
	
#From 4xn seqfreq matrix, calculate matrix of information contents
def seq_info_matrix(seq_freq_input, n):
	seqinfomatrix = np.zeros(seq_freq_input.shape)
	print("n = ", n)
	n = int(n)
	newFreq = float(round(1/(n+2), 5))
	print("n = ", n, "newFreq = ", newFreq)
	for x in range(np.shape(seq_freq_input)[0]):
		for y in range(np.shape(seq_freq_input)[1]):
			if seq_freq_input[x, y]==0:
				#newFreq = 0.002188
				info_content = 2 - (-1*math.log(newFreq, 2))
			else:
				#info_content = 2 + math.log(seq_freq_input[x, y], 2)
				info_content = 2 - (-1*math.log(seq_freq_input[x, y], 2))
			seqinfomatrix[x, y] = info_content
	return seqinfomatrix
	
#From 4xn matrix, calculate information content for alignment
def seq_info(sii):
	result = float(0)
	for colnbr in range(np.shape(sii)[1]):
		if np.sum(sii[:,colnbr]) == 0:
			result = result
		else:
			result = result + 2 - (-1*(
			sii[0, colnbr]*(0 if sii[0,colnbr]==0 else math.log(sii[0, colnbr], 2)) +
			sii[1, colnbr]*(0 if sii[1,colnbr]==0 else math.log(sii[1, colnbr], 2)) +
			sii[2, colnbr]*(0 if sii[2,colnbr]==0 else math.log(sii[2, colnbr], 2)) +
			sii[3, colnbr]*(0 if sii[3,colnbr]==0 else math.log(sii[3, colnbr], 2))))
	return result

def seq_info_sharp(sii):
	result = float(0)
	for colnbr in range(np.shape(sii)[1]):
		if np.sum(sii[:,colnbr]) == 0:
			result = result
		else:
			result = result + 2 - (-1*(
			(sii[0, colnbr]**3)*(0 if sii[0,colnbr]==0 else (math.log(sii[0, colnbr], 2)**4)) +
			(sii[1, colnbr]**3)*(0 if sii[1,colnbr]==0 else (math.log(sii[1, colnbr], 2)**4)) +
			(sii[2, colnbr]**3)*(0 if sii[2,colnbr]==0 else (math.log(sii[2, colnbr], 2)**4)) +
			(sii[3, colnbr]**3)*(0 if sii[3,colnbr]==0 else (math.log(sii[3, colnbr], 2)**4))))
	return result
	
def seq_info_sharper(sii):
	result = float(0)
	for colnbr in range(np.shape(sii)[1]):
		if np.sum(sii[:,colnbr]) == 0:
			result = result
		else:
			result = result + 2 - (-1*(
			(sii[0, colnbr]**6)*(0 if sii[0,colnbr]==0 else (math.log(sii[0, colnbr], 2)**8)) +
			(sii[1, colnbr]**6)*(0 if sii[1,colnbr]==0 else (math.log(sii[1, colnbr], 2)**8)) +
			(sii[2, colnbr]**6)*(0 if sii[2,colnbr]==0 else (math.log(sii[2, colnbr], 2)**8)) +
			(sii[3, colnbr]**6)*(0 if sii[3,colnbr]==0 else (math.log(sii[3, colnbr], 2)**8))))
	return result

def seq_info_enrich(sii):
	result = float(0)
	for colnbr in range(np.shape(sii)[1]):
		if np.sum(sii[:,colnbr]) == 0:
			result = result
		else:
			result = result + 2 - (-1*(
			sii[0, colnbr]*(0 if sii[0,colnbr]==0 else math.log(sii[0, colnbr], 2)) + ((sii[0, colnbr]/2)**2) +
			sii[1, colnbr]*(0 if sii[1,colnbr]==0 else math.log(sii[1, colnbr], 2)) + ((sii[1, colnbr]/2)**2) +
			sii[2, colnbr]*(0 if sii[2,colnbr]==0 else math.log(sii[2, colnbr], 2)) + ((sii[2, colnbr]/2)**2) +
			sii[3, colnbr]*(0 if sii[3,colnbr]==0 else math.log(sii[3, colnbr], 2)) + ((sii[3, colnbr]/2)**2)))
	return result
	
def compress_matrix(matrix_input):
#This will reduce the size of matrices to get rid of spaces from alignments
	gaps = np.zeros([matrix_input.shape[0], 2]) #stores register and total gaps
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
	output = np.zeros([matrix_input.shape[0], matrix_input.shape[1]-int(gaps[0,1])], dtype='str')
	for x in range(output.shape[0]):
		for y in range(output.shape[1]):
			output[x, y] = matrix_input[x, y+int(gaps[x,0])]
	return output
	
def malign_sharp(seq_matrix, seq_slide, seq_start, seq_stop, improvement_thresh, left_gap, seq_length):
	print("\nMalign: Sharp")
	start_time = time.clock()
	seq_matrix_w1 = seq_matrix[:,seq_start:seq_stop]
	seq_matrix_w2_size = [seq_matrix_w1.shape[0], seq_matrix_w1.shape[1]+2*seq_slide]
	seq_matrix_w2 = np.zeros(seq_matrix_w2_size, dtype='str')
	
	#Transferring sequences from seq.matrix.w1 to the middle of seq.matrix.w2
	for x in range(0, int(seq_matrix_w1.shape[0])):
		for y in range(0, int(seq_matrix_w1.shape[1])):
			seq_matrix_w2[x, y+seq_slide] = seq_matrix_w1[x,y]
					
	print(seq_matrix.shape[0], " sequences")
	#Heuristic start: pick five random sequences
	mat_heur = np.zeros([5, int(seq_matrix_w2_size[1])], dtype='str')
	for rowNum in range(mat_heur.shape[0]):
		randRowNum = random.randint(0, seq_matrix_w2.shape[1])
		mat_heur[rowNum,:]=seq_matrix_w2[randRowNum,:]
	#Heuristic start: brute force through all cases
	w_heur = np.zeros(mat_heur.shape, dtype='str')
	w_heur_best = np.zeros(mat_heur.shape, dtype='str')
	w_info_max = 0
	w_info_current = 0
	for s0 in range(2*seq_slide+1):
		for s1 in range(2*seq_slide+1):
			for s2 in range(2*seq_slide+1):
				for s3 in range(2*seq_slide+1):
					for s4 in range(2*seq_slide+1):
						for y in range(w_heur.shape[1]-2*seq_slide):
							w_heur[0, s0 + y] = mat_heur[0, y + seq_slide]
							w_heur[1, s1 + y] = mat_heur[1, y + seq_slide]
							w_heur[2, s2 + y] = mat_heur[2, y + seq_slide]
							w_heur[3, s3 + y] = mat_heur[3, y + seq_slide]
							w_heur[4, s4 + y] = mat_heur[4, y + seq_slide]
						w_info_current = seq_info_sharp(seq_freq(w_heur[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
						if w_info_current>w_info_max: #Pick each better alignment found
							w_info_max = w_info_current
							w_heur_best = w_heur
						if w_info_current==w_info_max and random.randint(0,1)==1:#randomly pick between two equal alignments found
							w_info_max = w_info_current
							w_heur_best = w_heur
						w_heur = np.zeros(mat_heur.shape, dtype='str')
	#Heuristic start done! w_heur now contains the best consensus found from 5 random sequences
	endtime_time = time.clock()
	print ("Heuristic start done. ", int(endtime_time-start_time), " seconds elapsed.")
	#First pass will append individual sequences to the heuristic consensus above and optimize
	w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str')
	for x in range(w_heur_best.shape[0]):
		for y in range(w_heur_best.shape[1]):
			w_firstp[x, y] = w_heur_best[x, y]
	w_firstp_info = 0
	w_firstp_info_max = 0
	wbest_firstp = np.zeros(w_firstp.shape, dtype='str')
	for s in range(2*seq_slide+1):
		w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str') #For some reason I have to reinitialize the matrix, and cannot say w_firstp[-1,:] = "". Check with Sumeet
		for x in range(w_heur_best.shape[0]):
			for y in range(w_heur_best.shape[1]):
				w_firstp[x, y] = w_heur_best[x, y]
		for y in range(w_firstp.shape[1]-2*seq_slide):
			w_firstp[-1, s + y] = seq_matrix_w2[0, y + seq_slide]
		w_firstp_info = seq_info_sharp(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))   ###
		if w_firstp_info>w_firstp_info_max:
			wbest_firstp = w_firstp
			w_firstp_info_max = w_firstp_info
			
	w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
	for x in range(wbest_firstp.shape[0]):
		for y in range(wbest_firstp.shape[1]):
			w_firstp[x, y] = wbest_firstp[x, y]
	for z in range(1,seq_matrix_w2.shape[0]):
		w_firstp_info_max = 0
		w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
		for s in range(2*seq_slide+1):
			w_firstp = np.zeros([w_firstp.shape[0], w_firstp.shape[1]], dtype='str')
			for x in range(wbest_firstp.shape[0]):
				for y in range(wbest_firstp.shape[1]):
					w_firstp[x, y] = wbest_firstp[x, y]
			for n in range(w_firstp.shape[1]):
				w_firstp[-1, n] = ""
			for n in range(w_firstp.shape[1]-2*seq_slide):
				w_firstp[-1, s + n] = seq_matrix_w1[z, n]
			w_firstp_info = seq_info_sharp(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
			if w_firstp_info>w_firstp_info_max:
				wbest_firstp = np.zeros(w_firstp.shape)
				wbest_firstp = w_firstp
				w_firstp_info_max = w_firstp_info
	endtime_time = time.clock()
	print("0th Pass Done. ", int(endtime_time-start_time), " seconds elapsed")
	
	#Remaining passes will shuffle individual sequences in preceeding passes
	w_nextp = np.zeros([wbest_firstp.shape[0]-5, wbest_firstp.shape[1]], dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = wbest_firstp[x+5, y]
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_1 = seq_info_sharp(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		pass_counter += 1
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""	
				for y in range(seq_matrix_w1.shape[1]):
					w_nextp[line, s + y] = seq_matrix_w1[line, y]
				line_info = seq_info_sharp(seq_freq(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
			wnextp_info_2 = seq_info_sharp(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign: Sharp is stopping.")
	return(w_nextp_best)

def malign_sharper(seq_matrix, seq_slide, seq_start, seq_stop, improvement_thresh, left_gap, seq_length):
	print("Malign: Sharper!!")
	start_time = time.clock()
	seq_matrix_w1 = seq_matrix[:,seq_start:seq_stop]
	seq_matrix_w2_size = [seq_matrix_w1.shape[0], seq_matrix_w1.shape[1]+2*seq_slide]
	seq_matrix_w2 = np.zeros(seq_matrix_w2_size, dtype='str')
	
	#Transferring sequences from seq.matrix.w1 to the middle of seq.matrix.w2
	for x in range(0, int(seq_matrix_w1.shape[0])):
		for y in range(0, int(seq_matrix_w1.shape[1])):
			seq_matrix_w2[x, y+seq_slide] = seq_matrix_w1[x,y]
					
	print(seq_matrix.shape[0], " sequences")
	#Heuristic start: pick five random sequences
	mat_heur = np.zeros([5, int(seq_matrix_w2_size[1])], dtype='str')
	for rowNum in range(mat_heur.shape[0]):
		randRowNum = random.randint(0, seq_matrix_w2.shape[1])
		mat_heur[rowNum,:]=seq_matrix_w2[randRowNum,:]
	#Heuristic start: brute force through all cases
	w_heur = np.zeros(mat_heur.shape, dtype='str')
	w_heur_best = np.zeros(mat_heur.shape, dtype='str')
	w_info_max = 0
	w_info_current = 0
	for s0 in range(2*seq_slide+1):
		for s1 in range(2*seq_slide+1):
			for s2 in range(2*seq_slide+1):
				for s3 in range(2*seq_slide+1):
					for s4 in range(2*seq_slide+1):
						for y in range(w_heur.shape[1]-2*seq_slide):
							w_heur[0, s0 + y] = mat_heur[0, y + seq_slide]
							w_heur[1, s1 + y] = mat_heur[1, y + seq_slide]
							w_heur[2, s2 + y] = mat_heur[2, y + seq_slide]
							w_heur[3, s3 + y] = mat_heur[3, y + seq_slide]
							w_heur[4, s4 + y] = mat_heur[4, y + seq_slide]
						w_info_current = seq_info_sharper(seq_freq(w_heur[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
						if w_info_current>w_info_max: #Pick each better alignment found
							w_info_max = w_info_current
							w_heur_best = w_heur
						if w_info_current==w_info_max and random.randint(0,1)==1:#randomly pick between two equal alignments found
							w_info_max = w_info_current
							w_heur_best = w_heur
						w_heur = np.zeros(mat_heur.shape, dtype='str')
	#Heuristic start done! w_heur now contains the best consensus found from 5 random sequences
	endtime_time = time.clock()
	print ("Heuristic start done. ", int(endtime_time-start_time), " seconds elapsed.")
	#First pass will append individual sequences to the heuristic consensus above and optimize
	w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str')
	for x in range(w_heur_best.shape[0]):
		for y in range(w_heur_best.shape[1]):
			w_firstp[x, y] = w_heur_best[x, y]
	w_firstp_info = 0
	w_firstp_info_max = 0
	wbest_firstp = np.zeros(w_firstp.shape, dtype='str')
	for s in range(2*seq_slide+1):
		w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str') #For some reason I have to reinitialize the matrix, and cannot say w_firstp[-1,:] = "". Check with Sumeet
		for x in range(w_heur_best.shape[0]):
			for y in range(w_heur_best.shape[1]):
				w_firstp[x, y] = w_heur_best[x, y]
		for y in range(w_firstp.shape[1]-2*seq_slide):
			w_firstp[-1, s + y] = seq_matrix_w2[0, y + seq_slide]
		w_firstp_info = seq_info_sharper(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))   ###
		if w_firstp_info>w_firstp_info_max:
			wbest_firstp = w_firstp
			w_firstp_info_max = w_firstp_info
			
	w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
	for x in range(wbest_firstp.shape[0]):
		for y in range(wbest_firstp.shape[1]):
			w_firstp[x, y] = wbest_firstp[x, y]
	for z in range(1,seq_matrix_w2.shape[0]):
		w_firstp_info_max = 0
		w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
		for s in range(2*seq_slide+1):
			w_firstp = np.zeros([w_firstp.shape[0], w_firstp.shape[1]], dtype='str')
			for x in range(wbest_firstp.shape[0]):
				for y in range(wbest_firstp.shape[1]):
					w_firstp[x, y] = wbest_firstp[x, y]
			for n in range(w_firstp.shape[1]):
				w_firstp[-1, n] = ""
			for n in range(w_firstp.shape[1]-2*seq_slide):
				w_firstp[-1, s + n] = seq_matrix_w1[z, n]
			w_firstp_info = seq_info_sharper(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
			if w_firstp_info>w_firstp_info_max:
				wbest_firstp = np.zeros(w_firstp.shape)
				wbest_firstp = w_firstp
				w_firstp_info_max = w_firstp_info
	endtime_time = time.clock()
	print("0th Pass Done. ", int(endtime_time-start_time), " seconds elapsed")
	
	#Remaining passes will shuffle individual sequences in preceeding passes
	w_nextp = np.zeros([wbest_firstp.shape[0]-5, wbest_firstp.shape[1]], dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = wbest_firstp[x+5, y]
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_1 = seq_info_sharper(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		pass_counter += 1
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""	
				for y in range(seq_matrix_w1.shape[1]):
					w_nextp[line, s + y] = seq_matrix_w1[line, y]
				line_info = seq_info_sharp(seq_freq(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
			wnextp_info_2 = seq_info_sharper(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign: SharpER is stopping.")
	return(w_nextp_best)
	
def malign(seq_matrix, seq_slide, seq_start, seq_stop, improvement_thresh, left_gap, seq_length):
	print("Malign")
	start_time = time.clock()
	seq_matrix_w1 = seq_matrix[:,seq_start:seq_stop]
	seq_matrix_w2_size = [seq_matrix_w1.shape[0], seq_matrix_w1.shape[1]+2*seq_slide]
	seq_matrix_w2 = np.zeros(seq_matrix_w2_size, dtype='str')
	
	#Transferring sequences from seq.matrix.w1 to the middle of seq.matrix.w2
	for x in range(0, int(seq_matrix_w1.shape[0])):
		for y in range(0, int(seq_matrix_w1.shape[1])):
			seq_matrix_w2[x, y+seq_slide] = seq_matrix_w1[x,y]
					
	print(seq_matrix.shape[0], " sequences")
	#Heuristic start: pick five random sequences
	mat_heur = np.zeros([5, int(seq_matrix_w2_size[1])], dtype='str')
	for rowNum in range(mat_heur.shape[0]):
		randRowNum = random.randint(0, seq_matrix_w2.shape[1])
		mat_heur[rowNum,:]=seq_matrix_w2[randRowNum,:]
	#Heuristic start: brute force through all cases
	w_heur = np.zeros(mat_heur.shape, dtype='str')
	w_heur_best = np.zeros(mat_heur.shape, dtype='str')
	w_info_max = 0
	w_info_current = 0
	for s0 in range(2*seq_slide+1):
		for s1 in range(2*seq_slide+1):
			for s2 in range(2*seq_slide+1):
				for s3 in range(2*seq_slide+1):
					for s4 in range(2*seq_slide+1):
						for y in range(w_heur.shape[1]-2*seq_slide):
							w_heur[0, s0 + y] = mat_heur[0, y + seq_slide]
							w_heur[1, s1 + y] = mat_heur[1, y + seq_slide]
							w_heur[2, s2 + y] = mat_heur[2, y + seq_slide]
							w_heur[3, s3 + y] = mat_heur[3, y + seq_slide]
							w_heur[4, s4 + y] = mat_heur[4, y + seq_slide]
						w_info_current = seq_info(seq_freq(w_heur[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
						if w_info_current>w_info_max: #Pick each better alignment found
							w_info_max = w_info_current
							w_heur_best = w_heur
						if w_info_current==w_info_max and random.randint(0,1)==1:#randomly pick between two equal alignments found
							w_info_max = w_info_current
							w_heur_best = w_heur
						w_heur = np.zeros(mat_heur.shape, dtype='str')
	#Heuristic start done! w_heur now contains the best consensus found from 5 random sequences
	endtime_time = time.clock()
	print ("Heuristic start done. ", int(endtime_time-start_time), " seconds elapsed.")
	#First pass will append individual sequences to the heuristic consensus above and optimize
	w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str')
	for x in range(w_heur_best.shape[0]):
		for y in range(w_heur_best.shape[1]):
			w_firstp[x, y] = w_heur_best[x, y]
	w_firstp_info = 0
	w_firstp_info_max = 0
	wbest_firstp = np.zeros(w_firstp.shape, dtype='str')
	for s in range(2*seq_slide+1):
		w_firstp = np.zeros([w_heur_best.shape[0]+1, w_heur_best.shape[1]], dtype='str') #For some reason I have to reinitialize the matrix, and cannot say w_firstp[-1,:] = "". Check with Sumeet
		for x in range(w_heur_best.shape[0]):
			for y in range(w_heur_best.shape[1]):
				w_firstp[x, y] = w_heur_best[x, y]
		for y in range(w_firstp.shape[1]-2*seq_slide):
			w_firstp[-1, s + y] = seq_matrix_w2[0, y + seq_slide]
		w_firstp_info = seq_info(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))   ###
		if w_firstp_info>w_firstp_info_max:
			wbest_firstp = w_firstp
			w_firstp_info_max = w_firstp_info
			
	w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
	for x in range(wbest_firstp.shape[0]):
		for y in range(wbest_firstp.shape[1]):
			w_firstp[x, y] = wbest_firstp[x, y]
	for z in range(1,seq_matrix_w2.shape[0]):
		w_firstp_info_max = 0
		w_firstp = np.zeros([wbest_firstp.shape[0]+1, wbest_firstp.shape[1]], dtype='str')
		for s in range(2*seq_slide+1):
			w_firstp = np.zeros([w_firstp.shape[0], w_firstp.shape[1]], dtype='str')
			for x in range(wbest_firstp.shape[0]):
				for y in range(wbest_firstp.shape[1]):
					w_firstp[x, y] = wbest_firstp[x, y]
			for n in range(w_firstp.shape[1]):
				w_firstp[-1, n] = ""
			for n in range(w_firstp.shape[1]-2*seq_slide):
				w_firstp[-1, s + n] = seq_matrix_w1[z, n]
			w_firstp_info = seq_info(seq_freq(w_firstp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
			if w_firstp_info>w_firstp_info_max:
				wbest_firstp = np.zeros(w_firstp.shape)
				wbest_firstp = w_firstp
				w_firstp_info_max = w_firstp_info
	endtime_time = time.clock()
	print("0th Pass Done. ", int(endtime_time-start_time), " seconds elapsed")
	
	#Remaining passes will shuffle individual sequences in preceeding passes
	w_nextp = np.zeros([wbest_firstp.shape[0]-5, wbest_firstp.shape[1]], dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = wbest_firstp[x+5, y]
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_1 = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		pass_counter += 1
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""	
				for y in range(seq_matrix_w1.shape[1]):
					w_nextp[line, s + y] = seq_matrix_w1[line, y]
				line_info = seq_info(seq_freq(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
			wnextp_info_2 = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign is stopping.")
	return(w_nextp_best)
	
def malign_iterate_sharp(matrix_input, improvement_thresh, seq_slide, left_gap, seq_length):
	print("Malign Iterative: Sharp")
	start_time = time.clock()
	gaps = np.zeros([matrix_input.shape[0], 2]) #stores register and total gaps
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
		#print(line, register, line_gaps)
	condense = np.zeros([matrix_input.shape[0], matrix_input.shape[1]-2*seq_slide], dtype='str')
	print(condense.shape)
	for x in range(condense.shape[0]):
		for y in range(condense.shape[1]):
			#print(x)
			#print(y)
			newY = y+int(gaps[x,0])
			#print(newY)
			condense[x, y] = matrix_input[x, y+int(gaps[x,0])]
	#seq_slide = gaps[0,1]
	
	print(matrix_input.shape[0], " sequences")
	w_nextp = np.zeros(matrix_input.shape, dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = matrix_input[x, y]
	wnextp_info_1 = seq_info_sharp(seq_freq_orig(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""			
				for y in range(condense.shape[1]):
					w_nextp[line, s + y] = condense[line, y]
				line_info = seq_info_sharp(seq_freq(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
		wnextp_info_2 = seq_info_sharp(seq_freq_orig(w_nextp_best[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		pass_counter += 1
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
		if wnextp_info_2 == wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign Iterative: Sharp is stopping.")
	return(w_nextp_best)
	
def malign_iterate_sharper(matrix_input, improvement_thresh, seq_slide, left_gap, seq_length):
	print("Malign Iterative: SharpER")
	print("matrix input shape = ", matrix_input.shape)
	print("seq_length = ", seq_length)
	with open("matrix_input.txt","w") as f:
		for x in range(matrix_input.shape[0]):
			for y in range(matrix_input.shape[1]):
				if matrix_input[x,y] =='':
					f.write("-")
				else:
					f.write(matrix_input[x,y])
			f.write("\n")
	start_time = time.clock()
	gaps = np.zeros([matrix_input.shape[0], 2]) #stores register and total gaps
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
	condense = np.zeros([matrix_input.shape[0], matrix_input.shape[1]-2*seq_slide], dtype='str')
	print("condense shape = ", condense.shape)
	for x in range(condense.shape[0]):
		### x = row number in condense
		#print("gap info for current line: ", gaps[x,0], gaps[x,1])
		for y in range(condense.shape[1]):
			#### y = col position in condense, this should be <= to 
			#### the sequence length without any gaps, i.e. the seq_length variable
			newY = y+int(gaps[x,0])
			#print(x,y,newY)
			condense[x, y] = matrix_input[x, y+int(gaps[x,0])]
	#seq_slide = gaps[0,1]
	
	print(matrix_input.shape[0], " sequences")
	w_nextp = np.zeros(matrix_input.shape, dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = matrix_input[x, y]
	wnextp_info_1 = seq_info_sharper(seq_freq_orig(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""			
				for y in range(condense.shape[1]):
					w_nextp[line, s + y] = condense[line, y]
				line_info = seq_info_sharper(seq_freq(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
		wnextp_info_2 = seq_info_sharper(seq_freq_orig(w_nextp_best[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		pass_counter += 1
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
		if wnextp_info_2 == wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign Iterative: SharpER is stopping.")
	return(w_nextp_best)
	
def malign_iterate(matrix_input, improvement_thresh, seq_slide, left_gap, seq_length):
	print("Malign Iterative")
	start_time = time.clock()
	gaps = np.zeros([matrix_input.shape[0], 2]) #stores register and total gaps
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
	condense = np.zeros([matrix_input.shape[0], matrix_input.shape[1]-2*seq_slide], dtype='str')
	for x in range(condense.shape[0]):
		for y in range(condense.shape[1]):
			condense[x, y] = matrix_input[x, y+int(gaps[x,0])]
	#seq_slide = gaps[0,1]
	
	print(matrix_input.shape[0], " sequences")
	w_nextp = np.zeros(matrix_input.shape, dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = matrix_input[x, y]
	wnextp_info_1 = seq_info(seq_freq_orig(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""			
				for y in range(condense.shape[1]):
					w_nextp[line, s + y] = condense[line, y]
				line_info = seq_info(seq_freq(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
		wnextp_info_2 = seq_info(seq_freq_orig(w_nextp_best[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		pass_counter += 1
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
		if wnextp_info_2 == wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign Iterative is stopping.")
	return(w_nextp_best)
	
def malign_iterate_enrich(matrix_input, improvement_thresh, seq_slide, left_gap, seq_length):
	print("Malign Iterative: Enriching")
	start_time = time.clock()
	gaps = np.zeros([matrix_input.shape[0], 2]) #stores register and total gaps
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
	condense = np.zeros([matrix_input.shape[0], matrix_input.shape[1]-2*seq_slide], dtype='str')
	for x in range(condense.shape[0]):
		for y in range(condense.shape[1]):
			condense[x, y] = matrix_input[x, y+int(gaps[x,0])]
	#seq_slide = gaps[0,1]
	
	print(matrix_input.shape[0], " sequences")
	w_nextp = np.zeros(matrix_input.shape, dtype='str')
	for x in range(w_nextp.shape[0]):
		for y in range(w_nextp.shape[1]):
			w_nextp[x, y] = matrix_input[x, y]
	wnextp_info_1 = seq_info_enrich(seq_freq_orig(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
	print_info = seq_info(seq_freq_orig(w_nextp[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
	print("Starting information content is: ", print_info)
	wnextp_info_2 = wnextp_info_1 #Initialize variable
	pass_counter = 0
	w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
	while True:
		wnextp_info_1 = wnextp_info_2 #Set old value to 'new' value from previous pass
		for line in range(w_nextp.shape[0]):
			line_info = 0
			line_info_max = 0
			w_nextp_best = np.zeros(w_nextp.shape, dtype='str')
			for s in range(2*seq_slide+1):
				for y in range(w_nextp.shape[1]):
					w_nextp[line, y] = ""			
				for y in range(condense.shape[1]):
					w_nextp[line, s + y] = condense[line, y]
				line_info = seq_info_enrich(seq_freq(w_nextp[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
				if line_info>line_info_max:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
				if line_info==line_info_max and random.randint(0,1)==1:
					line_info_max = line_info
					for v1 in range(w_nextp_best.shape[0]):
						for v2 in range(w_nextp_best.shape[1]):
							w_nextp_best[v1, v2] = w_nextp[v1, v2]
			w_nextp = w_nextp_best
		wnextp_info_2 = seq_info_enrich(seq_freq_orig(w_nextp_best[:, int(seq_slide+left_gap):int(seq_slide+left_gap+seq_length)])) ###
		print_info = seq_info(seq_freq_orig(w_nextp_best[:, seq_slide+left_gap:seq_slide+left_gap+seq_length]))
		pass_counter += 1
		endtime_time = time.clock()
		print("Pass # ", int(pass_counter), " information content is: ", print_info, ". ", int(endtime_time-start_time), " seconds elapsed")
		if wnextp_info_2<wnextp_info_1*improvement_thresh and wnextp_info_2>wnextp_info_1:
			break
		if wnextp_info_2 == wnextp_info_1:
			break
	print("Sequence improvement threshold reached. Malign Iterative: Enriching is stopping.")
	return(w_nextp_best)
	
def align_to_cons(aligned_input, unaligned_input, left_gap, seq_length):
	seq_slide = int((aligned_input.shape[1] - unaligned_input.shape[1])/2)
	scores = seq_info_matrix(seq_freq(aligned_input), aligned_input.shape[0])
	alignment_output = np.zeros([unaligned_input.shape[0], aligned_input.shape[1]], dtype='str')
	for line in range(unaligned_input.shape[0]):
		register_info_max = 0
		for register in range(2*seq_slide+1):
			register_info = 0
			line_seq_buffer = np.zeros([1, aligned_input.shape[1]], dtype='str')
			for y in range(unaligned_input.shape[1]):
				line_seq_buffer[0, register + y] = unaligned_input[line, y]
			#evaluate information of register
			for y in range(seq_slide+left_gap, seq_slide+left_gap+seq_length):
				score_pull = 0
				if line_seq_buffer[0, y]=='A':
					score_pull = scores[0, y]
				elif line_seq_buffer[0, y]=='T':
					score_pull = scores[1, y]
				elif line_seq_buffer[0, y]=='G':
					score_pull = scores[2, y]
				elif line_seq_buffer[0, y]=='C':
					score_pull = scores[3, y]
				else:
					score_pull = 0
				register_info = register_info + score_pull
			if register_info>register_info_max:
				register_info_max = register_info
				for y in range(line_seq_buffer.shape[1]):
					alignment_output[line, y] = line_seq_buffer[0, y]
	return alignment_output

def cherrypick_alignments(consensus_alignment, aligned_sequences, goodorbad, seq_slide, left_gap, seq_length, seqIDs):
#consensus alignment is the output from malign
#aligned_sequences is output from align_to_cons
#goodorbad is to be set to +1 for picking good sequences, -1 for bad sequences
	scores = seq_info_matrix(seq_freq_orig(consensus_alignment), aligned_sequences.shape[0])
	line_scores = np.zeros([aligned_sequences.shape[0], 1])
	for line in range(aligned_sequences.shape[0]):
		line_info = 0
		for y in range(seq_slide+left_gap, seq_slide+left_gap+seq_length):
			score_pull = 0
			if aligned_sequences[line, y]=='A':
				score_pull = scores[0, y]
			elif aligned_sequences[line, y]=='T':
				score_pull = scores[1, y]
			elif aligned_sequences[line, y]=='G':
				score_pull = scores[2, y]
			elif aligned_sequences[line, y]=='C':
				score_pull = scores[3, y]
			else:
				score_pull = 0
			line_info = line_info + score_pull
		line_scores[line, 0] = line_info
	output_good_alignments = np.zeros(aligned_sequences.shape, dtype='str')
	good_sequences = 0
	c = 0
	output_seqIDs = []
	output_lineScores = []
	for line in range(line_scores.shape[0]):
		if line_scores[line, 0]>0:
			output_seqIDs.append(seqIDs[c])
			output_lineScores.append(line_scores[line,0])
			for y in range(aligned_sequences.shape[1]):
				output_good_alignments[good_sequences, y] = aligned_sequences[line, y]
			good_sequences += 1
		c += 1
	output_matrix = np.zeros([good_sequences, output_good_alignments.shape[1]], dtype='str')
	for x in range(good_sequences):
		for y in range(output_good_alignments.shape[1]):
			output_matrix[x, y] = output_good_alignments[x, y]

	with open("EcoTSS.cherrypick.scores.txt","a") as f:
		f.write("\n".join(map(str, output_lineScores)))
		f.write("\n")
	with open("EcoTSS.cherrypick.seqIDs.txt","a") as f2:
		f2.write("\n".join(map(str, output_seqIDs)))
		f2.write("\n")
	return(output_matrix, output_seqIDs)

#Input alignment file from saved txt (use 'filename.ext' format)
def load_alignment(input):
	#Figure out the length of each line
	f = open(input, 'r')
	while True:
		temp = str(f.readline())
		if temp != '\n':
			break
	f.close()
	y = -1
	for chr in temp:
		y += 1
	
	#Figure out total number of rows in file
	f = open(input, 'r')
	x = 0
	while True:
		temp = str(f.readline())
		if temp != '\n' and temp != '':
			x += 1
		if temp == '':
			break
	f.close()
	
	array_1 = np.zeros([x, y], dtype = 'str')

	f = open(input, 'r')
	x1 = 0

	while True:
		temp = str(f.readline())
		if temp != '\n':
			y1 = 0
			for chr in temp:
				if chr != '\n':
					if chr == 'A' or chr == 'T' or chr == 'G' or chr == 'C':
						array_1[x1, y1] = chr
					y1 += 1
			x1 += 1
		if temp == '':
			break
	return array_1

#Save alignment file to saved txt
def save_alignment(array_alignment, filename):
	with open(filename,"w") as f:
		for x in range(array_alignment.shape[0]):
			for y in range(array_alignment.shape[1]):
				if str(array_alignment[x,y]) == '':
					f.write('-')
				else:
					f.write(array_alignment[x,y])
			f.write("\n")
			
def save_multiscanoutput(array_alignment, filename):
	with open(filename,"w") as f:
		for x in range(array_alignment.shape[0]):
			for y in range(array_alignment.shape[1]):
				f.write(str(array_alignment[x,y]))
				if y != array_alignment.shape[1]-1:
					f.write('\t')
			f.write("\n")
			
def load_multiscanoutput(filename):
	m = []
	with open(filename, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			m.append(list(row))
	m = np.array(m)
	m = m.astype(float)
	return m
			
def gap_list(matrix_input):
	gaps = np.zeros([matrix_input.shape[0], 2])
	for line in range(matrix_input.shape[0]):
		register = 0
		line_gaps = 0
		while True:
			if matrix_input[line, register]=='':
				register += 1
			else:
				break
		for y in range(matrix_input.shape[1]):
			if matrix_input[line, y]=='':
				line_gaps += 1
		gaps[line, 0] = register
		gaps[line, 1] = line_gaps
	return gaps
	
#Hex_10 and Hex_35 are both seq_freq matrices. seq_slide applies to -10 element search.
#spacing_set is ideally range(21, 27) (includes 21-26), and represents the spacings considered by the authors of Shultzaberger et al. I write the code such that it includes the max-spacing number
#This function returns the best -10 alignment hexamer
#Ideally pick the region including 59 - 96 from the promoter library(59 - 96)
#Outputs: -10 from left; -10 to -35; -10 Seq Info; -35 Seq Info; Gap Surprisal, Total info
#matrix_input, Hex_10, Hex_35, TSSloc: 100, sp3510min: -26, sp3510max: -21, sp10TSSmin: -14, sp10TSSmax: -8
def multiscan_promoter_10(matrix_input, Hex_10, Hex_35, TSSloc, sp3510min, sp3510max, sp10TSSmin, sp10TSSmax):
	scores_10 = seq_info_matrix(Hex_10, matrix_input.shape[0])
	scores_35 = seq_info_matrix(Hex_35, matrix_input.shape[0])
	output = np.zeros([matrix_input.shape[0], 6]) #Need to determine slots
	
	#Gap_Surprisal_Matrix_Denovo
	gap_surprisal_matrix = np.zeros([sp3510max-sp3510min+1, 2])
	for spacing in range(sp3510min, sp3510max+1):
		gap_surprisal_matrix[spacing-sp3510min, 0] = spacing
		gap_surprisal_matrix[spacing-sp3510min, 1] = 1+math.cos(2*math.pi*(spacing+23)/10.6)
	Totaln = 0
	for spacing in range(gap_surprisal_matrix.shape[0]):
		Totaln = Totaln + gap_surprisal_matrix[spacing, 1]
	for spacing in range(gap_surprisal_matrix.shape[0]):
		gap_surprisal_matrix[spacing, 1] = math.log(gap_surprisal_matrix[spacing, 1]/Totaln, 2)
	
	hexamer_test_10 = np.zeros([1, 6], dtype='str')
	hexamer_test_35 = np.zeros([1, 6], dtype='str')
	for line in range(matrix_input.shape[0]):
		info_total_best = -30
		for spacing_10_TSS in range(sp10TSSmin, sp10TSSmax+1):
			for y in range(hexamer_test_10.shape[1]):
				hexamer_test_10[0,y] = matrix_input[line, TSSloc+spacing_10_TSS-1+y]
			for spacing_35_10 in range(sp3510min, sp3510max+1):
				for y in range(hexamer_test_35.shape[1]):
					hexamer_test_35[0,y] = matrix_input[line, TSSloc+spacing_10_TSS+spacing_35_10-1+y]
				infogs = gap_surprisal_matrix[spacing_35_10-sp3510min, 1]
				info10total = 0
				info35total = 0
				for x in range(hexamer_test_10.shape[1]):
					score_pull = 0
					if hexamer_test_10[0, x]=='A':
						score_pull = scores_10[0, x]
					elif hexamer_test_10[0, x]=='T':
						score_pull = scores_10[1, x]
					elif hexamer_test_10[0, x]=='G':
						score_pull = scores_10[2, x]
					elif hexamer_test_10[0, x]=='C':
						score_pull = scores_10[3, x]
					else:
						score_pull = 0
					info10total = info10total + score_pull
				for y in range(hexamer_test_35.shape[1]):
					score_pull = 0
					if hexamer_test_35[0, y]=='A':
						score_pull = scores_35[0, y]
					elif hexamer_test_35[0, y]=='T':
						score_pull = scores_35[1, y]
					elif hexamer_test_35[0, y]=='G':
						score_pull = scores_35[2, y]
					elif hexamer_test_35[0, y]=='C':
						score_pull = scores_35[3, y]
					else:
						score_pull = 0
					info35total = info35total + score_pull
				info_total = infogs + info10total + info35total
				if info_total>info_total_best: #Update everything in the line
					info_total_best = info_total
					output[line, 0] = spacing_10_TSS #To be updated
					output[line, 1] = spacing_35_10
					output[line, 2] = info10total
					output[line, 3] = info35total
					output[line, 4] = infogs
					output[line, 5] = info_total
	return output
	
def multiscan_promoter_flex(matrix_input, Hex_10, Hex_35, TSSloc, sp3510min, sp3510max, gs3510, sp10TSSmin, sp10TSSmax, gs10TSS):
	scores_10 = seq_info_matrix(Hex_10, matrix_input.shape[0])
	scores_35 = seq_info_matrix(Hex_35, matrix_input.shape[0])
	output = np.zeros([matrix_input.shape[0], 6]) #Need to determine slots
	
	#Gap_Surprisal_Matrix_Imported_35_10
	gap_surprisal_matrix_3510 = np.zeros([sp3510max-sp3510min+1, 2])
	for spacing in range(gap_surprisal_matrix_3510.shape[0]):
		gap_surprisal_matrix_3510[spacing, 1] = gs3510[spacing]
		
	#Gap_Surprisal_Matrix_Imported_10_TSS
	gap_surprisal_matrix_10TSS = np.zeros([sp10TSSmax-sp10TSSmin+1, 2])
	for spacing in range(gap_surprisal_matrix_10TSS.shape[0]):
		gap_surprisal_matrix_10TSS[spacing, 1] = gs10TSS[spacing]
	
	hexamer_test_10 = np.zeros([1, 6], dtype='str')
	hexamer_test_35 = np.zeros([1, 6], dtype='str')
	for line in range(matrix_input.shape[0]):
		info_total_best = -30
		for spacing_10_TSS in range(sp10TSSmin, sp10TSSmax+1):
			for y in range(hexamer_test_10.shape[1]):
				hexamer_test_10[0,y] = matrix_input[line, TSSloc+spacing_10_TSS-1+y]
			for spacing_35_10 in range(sp3510min, sp3510max+1):
				for y in range(hexamer_test_35.shape[1]):
					hexamer_test_35[0,y] = matrix_input[line, TSSloc+spacing_10_TSS+spacing_35_10-1+y]
				infogs = gap_surprisal_matrix_3510[spacing_35_10-sp3510min, 1] + gap_surprisal_matrix_10TSS[spacing_10_TSS-sp10TSSmin, 1]
				info10total = 0
				info35total = 0
				for x in range(hexamer_test_10.shape[1]):
					score_pull = 0
					if hexamer_test_10[0, x]=='A':
						score_pull = scores_10[0, x]
					elif hexamer_test_10[0, x]=='T':
						score_pull = scores_10[1, x]
					elif hexamer_test_10[0, x]=='G':
						score_pull = scores_10[2, x]
					elif hexamer_test_10[0, x]=='C':
						score_pull = scores_10[3, x]
					else:
						score_pull = 0
					info10total = info10total + score_pull
				for y in range(hexamer_test_35.shape[1]):
					score_pull = 0
					if hexamer_test_35[0, y]=='A':
						score_pull = scores_35[0, y]
					elif hexamer_test_35[0, y]=='T':
						score_pull = scores_35[1, y]
					elif hexamer_test_35[0, y]=='G':
						score_pull = scores_35[2, y]
					elif hexamer_test_35[0, y]=='C':
						score_pull = scores_35[3, y]
					else:
						score_pull = 0
					info35total = info35total + score_pull
				info_total = infogs + info10total + info35total
				if info_total>info_total_best: #Update everything in the line
					info_total_best = info_total
					output[line, 0] = spacing_10_TSS #To be updated
					output[line, 1] = spacing_35_10
					output[line, 2] = info10total
					output[line, 3] = info35total
					output[line, 4] = infogs
					output[line, 5] = info_total
	return output
