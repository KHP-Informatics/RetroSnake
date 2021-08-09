# retroPlus

This is a SnakeMake pipeline based on RetroSeq https://github.com/tk2/RetroSeq,  a tool for discovery and genotyping of transposable element variants.
The pipeline runs RetroSeq, filters RetroSeq predictions, verifies the insertions by assembling the regions around each insertion and running the assembled contigs through RepeatMasker, and finally annotated the remaining high confidence insertions.


In this modular Snakemake pipeline, users can choose the level of filtering and verification of predicted insertions, as well as a few steps of downstream analysis: comparing the predictions with the known insertions to separate them into previously reported and novel insertions, as well as using AnnotSV to further annotate the insertions with genes and regulatory elements, and their potential clinical significance.


Snakemake pipelines are ideal for this scenario, as any result file can be created separately, and only steps necessary for the creation of the requested result file will be done. For example, if a user has already ran the pipeline up until ‘filterCalls’ – the step in which Retroseq predictions have been filtered to keep only the high quality ones, and now the pipeline is run again to request verified and annotated predictions, only the needed steps will be taken – in this case ‘verify’ and ‘annotate verified’. 

For example, if BAM file is paready present in the specified directory, the first step, CramtoBam conversion is skipped, and the BAM file goes directly to the retroseqDiscover and retroseqCall steps. 




As the snakemake pipeline is modular and providing different outputs, one way to call it is to specify the output file.  

For example,to produce a file with retroseq prediction which have been filtered and annotated, you would call it:

snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/LP6008463-DNA_G04.annotatedFiltered.tsv


For example,to produce a file with retroseq prediction which have been filtered and verified with the extra verification step (assembling the region and running throught RepeatMasker) and then annotated, you would call it:

snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/LP6008463-DNA_G04.annotatedVerified.tsv


To take filtered and verified predictions and to mark known and novel HERV insertions, you would call:
snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/ {LP6008463-DNA_G04.novelHitsFV.bed,LP6008463-DNA_G04.knownHitsFV.bed}

To run it on the cluster, you add:
 
--cluster "sbatch -p YOUR_PARTITION --mem-per-cpu=7G"

It can be run with multiple files simultaneously, snakemake will split it in different cores. For example, to get annotated and verified predictions for several files (BAM or CRAM), you would call it:

snakemake --use-conda --use-envmodules --cores 5 <RESULTS_DIRECTORY>/{BAM1_prefix,BAM2_prefix, BAM3_prefix}.annotatedVerified.tsv

To call in on the cluster

nohup snakemake -s SnakefileRetroPlus --use-conda --use-envmodules --cores 5 --cluster "sbatch -p brc --mem-per-cpu=7G" /mnt/lustre/groups/herv_project/snakemake/retroseq/results/{LP6008463-DNA_G04.annotatedFiltered.tsv,LP6008463-DNA_G04.annotatedFiltered.html,LP6008463-DNA_G04.annotatedVerified.tsv,LP6008463-DNA_G04.annotatedVerified.html,LP6008463-DNA_G04.novelHitsF.bed,LP6008463-DNA_G04.novelHitsFV.bed,LP6008463-DNA_G04.knownHitsF.bed,LP6008463-DNA_G04.knownHitsFV.bed}


# Installing Dependencies

## RepeatMasker
https://github.com/rmhubley/RepeatMasker


## AnnotSV 
Follow the instructions to download and install AnnotSV and knotAnnotSV.  This is only needed if you are going to run the annotation step of RetroPlus pipeline.  

https://lbgi.fr/AnnotSV/

https://github.com/mobidic/knotAnnotSV

Important:
After the installation, update the path of the two install directories (to AnnotSV and knotAnnotSV) in the config.yaml file.
