# local-alignment-for-protein-alignment
local alignment implementation with linear gap costs for protein alignment in python


################################################################################################################################

import sys

def blossum_row(aa, row):
  values = row [1:len(row)]
  blossum_row = {}
  for i in range(0,len(values)):
    blossum_row[aa[i]] = int(values[i])
  return blossum_row

def calc_score(blossum_matrix,matrix,i,j):
  diag_left = matrix[i-1][j-1] + blossum_matrix[matrix[i][0]][matrix[0][j]]
  left = matrix[i][j-1] - 10
  top = matrix[i-1][j] - 10
  max_value = max(0, diag_left, left, top)
  values = [0, diag_left, left, top]
  source = ['-', (i-1,j-1), (i,j-1), (i-1,j)]
  return (max_value, source[values.index(max_value)])
 

def find_best_value(matrix):
  best_value = 0
  i_best = 0
  j_best = 0
  for i in range(2,len(matrix)):
    for j in range(2,len(matrix[0])):
      if (matrix[i][j] > best_value):
        best_value = matrix[i][j]
        i_best = i
        j_best = j
  return (best_value, i_best,j_best)

def compute_local_alignment(seq_1, seq_2):
  blossum_matrix = {}
  file = open("BLOSUM62.txt", 'r')
  lines = file.readlines()
  aa = ""
  for i in range(0,len(lines)):
    if (lines[i][0] == '#'):
      continue
    elif (lines[i][0] == ' '):
      aa = "".join(lines[i].split())
    else:
      row = lines[i].split()
      blossum_r = blossum_row(aa,row)
      blossum_matrix[row[0]] = blossum_r
  
  ''' init matrix '''
  
  matrix = [[0 for x in range(len(seq_1)+2)] for x in range(len(seq_2)+2)]
  matrix[0][0] = ''
  matrix[0][1] = ''
  matrix[1][0] = ''
  for i in range(2,len(seq_2)+2):
    matrix[i][0] = seq_2[i-2] 
  for i in range(2,len(seq_1)+2):
    matrix[0][i] = seq_1[i-2]
  
  ''' init backtrace matrix ''''
  
  backtrace = [['-' for x in range(len(seq_1)+2)] for x in range(len(seq_2)+2)]
  backtrace[0][0] = ''
  backtrace[0][1] = ''
  backtrace[1][0] = ''
  for i in range(2,len(seq_2)+2):
    backtrace[i][0] = seq_2[i-2] 
  
  for i in range(2,len(seq_1)+2):
    backtrace[0][i] = seq_1[i-2] 
  
  for i in range(2,len(matrix)): 
    for j in range(2,len(matrix[0])):
      score, source = calc_score(blossum_matrix,matrix,i,j)
      matrix[i][j] = score
      backtrace [i][j] = source
  
  ''' backtrace'''
  
  a_seq_1 = ""
  a_seq_2 = ""
  best_value, i, j = find_best_value(matrix)
  
  while(backtrace[i][j] != '-'):
    i_new, j_new = backtrace[i][j]
    if (i-1 == i_new and j-1 == j_new):
      a_seq_1 = matrix[0][j] + a_seq_1
      a_seq_2 = matrix[i][0] + a_seq_2
      i = i_new
      j = j_new
    elif (i == i_new and j-1 == j_new):
      a_seq_1 = matrix[0][j] + a_seq_1
      a_seq_2 = "-" + a_seq_2
      i = i_new
      j = j_new
    elif (i-1==i_new and j==j_new):
      a_seq_1 = '-' + a_seq_1
      a_seq_2 = backtrace[i][0] + a_seq_2
      i = i_new
      j = j_new
   print result
  print(best_value)
  print(a_seq_1)
  print(a_seq_2)
      
''' main'''


seq_1 = sys.argv[1]
seq_2 = sys.argv[2]

compute_local_alignment(seq_1,seq_2)
