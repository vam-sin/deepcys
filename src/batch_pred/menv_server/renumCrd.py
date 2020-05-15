import math
import sys
input_file = sys.argv[1]
file_name = open(input_file,"r")
output_file = open("Crdresults.txt","w")
for line in file_name:
       output_file.write(line[6:11]+' '+line[13:17]+' '+line[62:70]+ '\n')
