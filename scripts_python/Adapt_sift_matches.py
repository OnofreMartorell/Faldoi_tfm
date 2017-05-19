print('Hello world')
f = open('../data/alley_2/matches.txt', 'r')
g = open('../data/alley_2/matches_sift.txt', 'w')

#file = f.read()

for line in f:
        values = line.split()
	values_new = [values[0], values[1],values[4],values[5]]
	str1 = ' '.join(str(e) for e in values_new)
	print(str1)
	g.write(str1 + '\n')

g.close()
