from sklearn import metrics

a = open("mod.txt")
b = open("actual.txt")

# Predictions from a (mod.txt)
pred = []
for i in a:
	i = i.replace("\n", "")
	# print(i)
	pred.append(i)

# Actual values from b (actual.txt)
act = []
j = 0
for i in b:
	i = i.replace("\n", "")
	# print(i)
	act.append(i)
	j += 1
	if j == 1: # Number predicted
		break

# print(act)

results = metrics.accuracy_score(act, pred)

print(results)
	
