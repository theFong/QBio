
fname = 'match.result.txt'
f = open(fname,'r')
lines = f.readlines()[4:]
str1 = ''
str2 = ''
check = 0
for line in lines:
	if check % 4 == 0:
		str1 += line

	if check % 4 == 2:
		str2 += line
	check += 1

f.close()
open('matchSeq.txt','w').write(str1)

fname = 'sampleoutput/Ebola_Zaire_vs_Reston_Comparison.result'
f = open(fname,'r')
lines = f.readlines()[4:]
str1 = ''
str2 = ''
check = 0
for line in lines:
	if check % 4 == 0:
		str1 += line

	if check % 4 == 2:
		str2 += line
	check += 1

f.close()
open('ebolaCorrect.txt','w').write(str1)