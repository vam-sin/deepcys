import math
import sys

B = []
for index, line in enumerate(open('uniqueMenvresidues.txt')):
	B.append(line.split()[1])

A = []
with open("uniqueresidues.txt") as f:
	for line in f:
		d={}
		key= line.split()[0]
		val= line.split()[1]
		d[int(key)] = val
		A.append(d)
f.close()

for index, item in enumerate(B):
	if B[index] == 'HSP':
		B[index] = 'HIS'

l =  len(A)
i = 0 
j = 1
C = []

while i < l:
	try :
		if(A[i][j] == B[i]):
			i+=1
			j+=1
			#print i, "success"
		else:
			#print A.pop(i)
			C.append(A.pop(i))
			j+=1
	except IndexError:
		k = 0
		for i in A:
			if k>= len(B):
				A.pop(k)
			k+=1
		break

#print C


all_keys = set().union(*(d.keys() for d in C))
Allkeys = list(all_keys)

#print Allkeys

wanted = Allkeys
output_file = open('outfile.txt','w')
with open('uniqueresidues.txt') as f:
	for line in f:
		if int(line[0:6]) not in wanted:
			output_file.write(line)
		else:
			print(line[7:])
