import math
import sys
input_file = sys.argv[1]
file_name = open(input_file,"r")
output_file = open("pythonresults.txt","w")
for line in file_name:
    if(line[0:4]=='ATOM'): 
       output_file.write(line[17:20]+' '+line[20:22]+' '+line[23:31]+ '\n')

