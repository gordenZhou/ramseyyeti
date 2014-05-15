ontoc = {}
r = 42
fp = open('../r55_42.txt')
for line in fp:
	line = line.rstrip('\n')
	if (not line.startswith('0') or line.startswith('1')):
		continue
	on = 0
	for i in range(0,r):
		if (line[i]=="1"):
			on += 1
	if not ontoc.has_key(on):
		ontoc[on] = 0
	ontoc[on] += 1
fp.close()

print ontoc