import sys

allContigs=[]
contigs=[]

inFile = open(sys.argv[1],'r')
element = sys.argv[2]
verificationLevel=sys.argv[3]
#print ("Checking the file " + sys.argv[1])
#inFile=open("chr4_9602976-9603226.fa.cap.contigs.out",'r')
#inFile=open("chr2_65138685-65138935.fa.cap.contigs.out",'r')

for line in inFile:
    bits=line.split()
    if len(bits)>0:
        if not (bits[4] in allContigs) and not (bits[4]=="query") and not (bits[4]=="sequence"):
            allContigs.append(bits[4])
        if line.find (element)>0:
             #32    6.1  0.0  0.0  Contig1       1    49  (491) C SVA_E     Retroposon/SVA  (1263)    119     71   1
            if float(bits[1])<10 and float(bits[2])<3 and float(bits[3])<3:
                contigs.append(bits[4])

inFile.close()
contigs.sort()
allContigs.sort()


#rules: either more than one contig, or if it is on, towards the top of the list
#/mnt/lustre/groups/herv_project/herv_pipeline_out/sampletemp2/LP6008118-DNA_D05/repeatmasker/chr9:132205151-132206151.fasta.cap.contigs.out
complete = sys.argv[1]
ind1=complete.find('/chr')
ind2=complete.find('.fasta.cap')
complete=complete[ind1+1:ind2]
chro=complete.split(':')[0]
start=complete.split(':')[1].split('-')[0]
stop=complete.split(':')[1].split('-')[1]

if verificationLevel == "high":

 if len (contigs)>1:
    #print (sys.argv[1])
    print (chro + "\t" + start + "\t" + stop) 
    #print ("hit")
    exit()
 if len(contigs)==1:
    position=allContigs.index(contigs[0])+1
    if position==1:
      #print (sys.argv[1])
      print (chro + "\t" + start + "\t" + stop) 
      exit()       
    if position/len(contigs) < 0.51:
        print (chro + "\t" + start + "\t" + stop) 
        #print (sys.argv[1])
        exit() 

else: #low verification level, if the element is found in ANY contig
   if len (contigs)>0:
      print (chro + "\t" + start + "\t" + stop)
      exit()
