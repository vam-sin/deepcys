f = open("data.txt", "r")
g = open("bench.txt", "w")

for i in f:
	i = i.split('\t')
	print(i)
	g.write(i[0] + ',' + i[1] + ',' + i[2].replace("\n", ""))
	g.write("\n")