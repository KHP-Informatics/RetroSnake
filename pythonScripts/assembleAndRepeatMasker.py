import os
import sys


inFile = open(sys.argv[1],"r")

print ("Writing to file " + sys.argv[2])
bamFile = sys.argv[3]
outDir = sys.argv[4]
scriptsDir=sys.argv[6]
element=sys.argv[7]
repeatMaskerPath = sys.argv[5]
for line in inFile:
    line = line.rstrip()
##    print (line)
    os.system("samtools view -b %s -o  %s/%s.bam %s" %(bamFile,outDir,line,line))
    os.system("samtools bam2fq %s/%s.bam > %s/%s.fq" %(outDir,line,outDir,line))
    os.system("python %s/FastqToFastaAndQual.py %s/%s" %(scriptsDir,outDir,line))
    os.system("cap3 %s/%s.fasta  > %s/%s.log" %(outDir,line,outDir,line))
    os.system("perl %sRepeatMasker %s/%s.fasta.cap.contigs -pa 2 -dir %s/repeatmasker" %(repeatMaskerPath,outDir,line, outDir)) 
    os.system("python %s/testForHervPresence.py %s/repeatmasker/%s.fasta.cap.contigs.out %s >> %s" %(scriptsDir,outDir,line,element,sys.argv[2]))

inFile.close()
