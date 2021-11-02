
import sys

inFile = open(sys.argv[1],'r')
outFile = open(sys.argv[2],'w')
outFile2 = open(sys.argv[3],'w')
#chr21	18047285	.	T	<INS:ME>	46	.	MEINFO=HERVK,18047285,18047286,NA;NOT_VALIDATED;SVTYPE=INS	GT:GQ:FL:SP:CLIP5:CLIP3	1/1:46:8:0:2:1
             
for line in inFile:             
    if line[0] == '#':
        continue
    
    ol = line
    line = line.rstrip()
    line = line.split()
    g = line[-1]
    fl = g.split(':')[2]
    fl = int(fl)
    gq = g.split(':')[1]
    gq = int(gq)

    #if ((fl == 6 and gq > 28) or (fl == 7 and gq > 20) or (fl == 8 and gq > 20)):
    if ((fl == 6 and gq > 28) or (fl == 7 and gq > 20) or (fl == 8 and gq > 10)):
            #outFile.write(ol)
            chro=line[0]
            start=line[1]
            end=int(line[1])+1000
            outFile.write(chro + ":" + start + "-" + str(end) + "\n")
            outFile2.write(chro + "\t" + start + "\t" + str(int(start)+1) + "\n")
    
inFile.close()
outFile.close()
outFile2.close()

