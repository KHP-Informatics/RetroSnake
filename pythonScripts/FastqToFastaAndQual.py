from Bio import SeqIO
import sys

#inputFastq=sys.argv[1]
prefix=sys.argv[1]
#SeqIO.convert(inputFastq, "fastq", "example.fasta", "fasta")
#SeqIO.convert(inputFastq, "fastq", "example.qual", "qual")
SeqIO.convert(prefix+".fq", "fastq", prefix+".fasta", "fasta")
SeqIO.convert(prefix+".fq", "fastq", prefix+".qual", "qual")
