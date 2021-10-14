# RetroSnake

This is a SnakeMake pipeline based on RetroSeq https://github.com/tk2/RetroSeq,  a tool for discovery of transposable element variants.
The pipeline runs RetroSeq configured to seach for HERV-K insertions (can look for any other mobile elements, see below how to change that), filters RetroSeq predictions, verifies the insertions by assembling the regions around each insertion and running the assembled contigs through RepeatMasker, and finally annotates insertions.

<img width="794" alt="image" src="https://user-images.githubusercontent.com/51742526/128677961-017e99c9-bdd0-4951-bbc3-e0057a3e1fe1.png">


Users can choose the level of filtering and verification of predicted insertions, as well as a few steps of downstream analysis: comparing the predictions with the known HERV-K insertions to separate them into previously reported and novel insertions, as well as using AnnotSV to further annotate the insertions with genes and regulatory elements, and their potential clinical significance.


As the snakemake pipeline is modular and providing different outputs, one way to call it is to specify the output file.

## Example 1

To produce a file with retroseq prediction which have been filtered and annotated, you would call it:

```
snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/<YOUR_SAMPLE_PREFIX>.annotatedFiltered.tsv
```
This would run retroseq, filter the results, and run AnnotSV to annotate them. The prerequisites are that either a <YOUR_SAMPLE_PREFIX>.bam or <YOUR_SAMPLE_PREFIX>.cram exist in the respective bam or cram directories specified in file config.yaml.


## Example 2

To produce a file with retroseq prediction which have been filtered and verified with the extra verification step (assembling the region and running throught RepeatMasker) and then annotated, you would call it:
```
snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/<YOUR_SAMPLE_PREFIX>.annotatedVerified.tsv
```

## Example 3

To take filtered and verified predictions and to mark known and novel HERV-K insertions, you would call:
```
snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/ {<YOUR_SAMPLE_PREFIX>.novelHitsFV.bed,<YOUR_SAMPLE_PREFIX>.knownHitsFV.bed}
```

## Example 4
To run it on the cluster, you add to the command line:
``` 
--cluster "sbatch -p YOUR_PARTITION --mem-per-cpu=3G"
```
## Example 5

It can be run with multiple files simultaneously, snakemake will split it in different cores. For example, to get annotated and verified predictions for several files (BAM or CRAM), you would call it:
```
snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/{BAM1_prefix,BAM2_prefix, BAM3_prefix}.annotatedVerified.tsv
```
## Example 6

To call in on the cluster
```
nohup snakemake --use-conda --use-envmodules --cores 5 --cluster "sbatch -p brc --mem-per-cpu=3G" /MY_OUTPUT_DIRECTORY/results/{BAM1_prefix.annotatedFiltered.tsv,BAM1_prefix.annotatedFiltered.html,BAM1_prefix.annotatedVerified.tsv,BAM1_prefix.annotatedVerified.html,BAM1_prefix.novelHitsF.bed,BAM1_prefix.novelHitsFV.bed,BAM1_prefix.knownHitsF.bed,BAM1_prefix.knownHitsFV.bed}
```
## Example 7

Alternatively, you can specify in the furst rule, all, the requested target files, and define the samples to run it on:

```
SAMPLES=["<YOUR_SAMPLE_PREFIX>"]

rule all:
   input:
        expand(outPath + "filter/{sample}.bed", sample=SAMPLES),
        expand(outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed",sample=SAMPLES),
        expand(outPath + "results/{sample}.knownHitsF.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.knownHitsFV.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.novelHitsF.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.novelHitsFV.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.annotatedFiltered.tsv", sample=SAMPLES),
        expand(outPath + "results/{sample}.annotatedFiltered.html", sample=SAMPLES),  
        expand(outPath + "results/{sample}.annotatedVerified.tsv", sample=SAMPLES),
        expand(outPath + "results/{sample}.annotatedVerified.html", sample=SAMPLES)
```
The snakemake can then be invoked without specifying the target files, for example:
```
nohup snakemake --use-conda --use-envmodules --cores 5 
```
or on the cluster:
```
nohup snakemake --use-conda --use-envmodules --cores 5 --cluster "sbatch -p MY_PARTITION --mem-per-cpu=3G" 
```
# Installing Dependencies

## Reference genome
If you already have a reference genome (and it is indexed), update the path to it in config.yaml
If you do not, run the script downloadHG19.sh in order to install it and index it.
```
bash downloadHG19.sh PATH_TO_DIR_TO_PLACE_REFERENCE_GENOME
```

## RepeatMasker
https://github.com/rmhubley/RepeatMasker
RepeatMasker is only needed for the verification step.  Without it you can still run RetroSeq, Filter insertions, mark known and novel insertions and run functional annotation. If you do not want the extra verification step, you do not need to install RepeatMasker.

Important:
After the installation of RepeatMasker, update the installation path in the config.yaml file.

## AnnotSV 
Follow the instructions to download and install AnnotSV and knotAnnotSV.  This is only needed if you are going to run the annotation step of RetroPlus pipeline.  

https://lbgi.fr/AnnotSV/

https://github.com/mobidic/knotAnnotSV

Important:
After the installation, update the path of the two install directories (to AnnotSV and knotAnnotSV) in the config.yaml file.

Unzip the file resources/insME_WithAF.bed.zip and place it in directory:
<ANNOT_SV>/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt

# Modifying the pipeline
## Change the transposable element to look for - the default is HERV-K

1.In config.yalm

HERVK_eref: /mnt/lustre/groups/herv_project/herv_pipeline_startfiles/hervk_eref.tab

hervk_eref.tab has the following contents:
HERVK	/mnt/lustre/groups/herv_project/herv_pipeline_startfiles/HERVK_ref.fasta

The first column is the name of transposable elementm the second its sequence in fasta file.

Replace the content of hervk_eref.tab file to point to FASTA of another transposable element. The name 'HERVK_eref' and 'hervk_eref.tab' can stay the same.

2.In config yaml
knownNR: /mnt/lustre/groups/herv_project/herv_pipeline_startfiles/list_of_known_integration_sitesNR.bed

Replace the file list_of_known_integration_sitesNR.bed with the list of known integration sites for your own transposable element. Thus is only needed for stepos X and Y.

3.In config yaml
element - replace with the name of TE, as defined in RepeatMasker libraries



## Retroseq filtering step - definition and adjustment

The filtering is done in script filterHighQualRetroseqForDownstream.py, found in folder pythonScripts.  It takes into account flags fl and gq. 

The actual filtering line is 
```
if ((fl == 6 and gq > 28) or (fl == 7 and gq > 20) or (fl == 8 and gq > 10)):
``` 
 To change the criteria, adjust the combination of flags to qualify an insertion in that line 
 
 
 # RepeatMasker output filtering
 
 This is achieved in script testForHervPresence.py.
 
 The first line of filtering - is a repetitive element found in a particular contig.
 
 
 Filtering criteria: 

% substitutions in matching region compared to the consensus < 10

% of bases opposite a gap in the query sequence (deleted bp) < 3

% of bases opposite a gap in the repeat consensus (inserted bp) < 3

These can be adjusted in line:
```
if float(bits[1])<10 and float(bits[2])<3 and float(bits[3])<3:
 ```           


The second is the verification level.  This can be changed by changing the parameter verificationLevel in rule verify.

If the verification level is "High" the rule is as follows:
If a TE (by default LTR5_Hs) is found either in more than one contig, OR if it is found only in one contig, that one would have to be in the top half of contigs as ranked by cap3;

If the verification level is "Low" (default), if there is a LTR5_Hs anywhere in any contig, it will pass the filter.
 
