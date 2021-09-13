import os
import sys


inFile = open(sys.argv[1],"r")

bamFile = sys.argv[2]
parts= bamFile.split('/')
bam=parts[-1][0:(len(parts[-1])-4)]
outDir = sys.argv[3]
outDirRep = outDir + "/" + bam
os.system("mkdir -p %s" %(outDirRep))
scriptsDir=sys.argv[5]
element=sys.argv[6]
repeatMaskerPath = sys.argv[4]
verificationLevel=sys.argv[7]
outFile=sys.argv[8]
for line in inFile:
    line = line.rstrip()
    os.system("samtools view -b %s -o  %s/%s.bam %s" %(bamFile,outDirRep,line,line))
    os.system("samtools bam2fq %s/%s.bam > %s/%s.fq" %(outDirRep,line,outDirRep,line))
    os.system("python %s/FastqToFastaAndQual.py %s/%s" %(scriptsDir,outDirRep,line))
    os.system("cap3 %s/%s.fasta  > %s/%s.log" %(outDirRep,line,outDirRep,line))
    os.system("perl %sRepeatMasker %s/%s.fasta.cap.contigs -pa 2 -dir %s/repeatmasker" %(repeatMaskerPath,outDirRep,line, outDirRep)) 
    os.system("python %s/testForHervPresence.py %s/repeatmasker/%s.fasta.cap.contigs.out %s %s >> %s" %(scriptsDir,outDirRep,line,element,verificationLevel,outFile))

inFile.close()

