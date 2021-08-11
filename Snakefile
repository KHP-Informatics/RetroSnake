configfile: "config.yaml"

SAMPLES = ["BAM1_Prefix","BAM2_Prefix","BAM3_Prefix"]
outPath = config["outPath"]
bamPath = config["bamPath"]
cramPath = config["cramPath"]


rule all:
   input:
        expand(outPath + "results/{sample}.knownHitsF.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.knownHitsFV.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.novelHitsF.bed", sample=SAMPLES),
        expand(outPath + "results/{sample}.novelHitsFV.bed", sample=SAMPLES)



rule CramToBam:
    input:
        cram_file=cramPath + "{sample}.cram"
    output:
        temp(bamPath + "{sample}.bam"),
        temp(bamPath + "{sample}.bam.bai")
    benchmark:
        "benchmarks/{sample}.CramToBam.benchmark.tsv"
    conda:
        "envs/sam_only.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        samtools view -b -h -@ {threads} -T {config[refHg19]} -o {output[0]} {input.cram_file}
        samtools index -@ {threads} {output[0]}
        """

rule retroseqDiscover:
    input:
        bamPath + "{sample}.bam",
        bamPath + "{sample}.bam.bai"
    output:
        outPath + "discover/{sample}.bed"
    threads: 8
    log:
        "logs/discover/{sample}.log"
    
    params:
        identity=80
    benchmark:
       #repeat("benchmarks/{sample}.retroseqDiscover.benchmark.txt",3)
       "benchmarks/{sample}.retroseqDiscover.benchmark.txt"
    conda:
       "envs/retroseq.yaml"
    shell:
        "perl {config[retroPath]}/retroseq.pl -discover -bam {input[0]} -output {output} -eref {config[HERVK_eref]} -id {params.identity}"
 
rule retroseqCall:
    input:
        bam=bamPath + "{sample}.bam",
        bai=bamPath + "{sample}.bam.bai",
        discover=outPath + "discover/{sample}.bed"
    output:
        outPath + "call/{sample}.vcf"
    threads: 8
    benchmark:
       "benchmarks/{sample}.retroCall.benchmark.txt"
    conda:
       "envs/retroseq.yaml"
    log:
        "logs/call/{sample}.log"
    shell:
        "perl {config[retroPath]}/retroseq.pl -call -bam {input.bam} -input {input.discover} -ref {config[refHg19]} -output {output}"

rule filterCalls:
   input:
      outPath + "call/{sample}.vcf"
   output:
      outPath + "filter/{sample}.pos",
      outPath + "filter/{sample}.bed"
   log:
        "logs/filter/{sample}.log"
   benchmark:
       "benchmarks/{sample}.filterCalls.benchmark.txt"
   shell:
       "python {config[pythonScripts]}/filterHighQualRetroseqForDownstream.py  {input} {output}"


rule verify:
    input:
      outPath + "filter/{sample}.pos",
      bamPath + "{sample}.bam",
      bamPath + "{sample}.bam.bai"
    output:
      outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed"
    benchmark:
      "benchmarks/{sample}.verify.tsv"
    conda:
       "envs/samtools.yaml"
    log:
        "logs/call/{sample}.log"
    shell:
      """
      python {config[pythonScripts]}/assembleAndRepeatMasker.py {input[0]} {output} {config[bamPath]}{wildcards.sample}.bam {config[outPath]} {config[RepeatMaskerPath]} {config[pythonScripts]} {config[element]} 
      """ 

rule markKnownFiltered:
    input:
      outPath + "filter/{sample}.bed"
    output:
       outPath + "results/{sample}.knownHitsF.bed"
    conda:
      "envs/bedtools.yaml"
    benchmark:
       "benchmarks/{sample}.markKnownFiltered.benchmark.txt"
    shell:
       "bedtools window -w 100 -c -a {config[knownNR]} -b  {input}  > {output}"

rule markNovelFiltered:
    input:
       outPath + "filter/{sample}.bed"
    output:
       outPath + "results/{sample}.novelHitsF.bed"
    conda:
      "envs/bedtools.yaml"
    benchmark:
       "benchmarks/{sample}.markNovelFiltered.benchmark.txt"
    shell:
        """
        bedtools sort -i  {input} >  {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        bedtools window -w 2000 -v -a {config[outPath]}novel/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output}  
        rm {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        """

rule markKnownVerified:
    input:
      outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed"
    output:
       outPath + "results/{sample}.knownHitsFV.bed"
    conda:
      "envs/bedtools.yaml"
    benchmark:
       "benchmarks/{sample}.markKnownVerified.benchmark.txt"
    shell:
       "bedtools window -w 100 -c -a {config[knownNR]} -b  {input}  > {output}"

rule markNovelVerified:
    input:
       outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed"
    output:
       outPath + "results/{sample}.novelHitsFV.bed"
    conda:
      "envs/bedtools.yaml"
    benchmark:
       "benchmarks/{sample}.markNovelVerified.benchmark.txt"
    shell:
        """
        bedtools sort -i  {input} >  {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        bedtools window -w 2000 -v -a {config[outPath]}novel/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output}  
        rm {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        """

rule annotateFiltered:
     input:
       outPath + "filter/{sample}.bed" 
     output:
       outPath + "results/{sample}.annotatedFiltered.tsv",
       outPath + "results/{sample}.annotatedFiltered.html" 
     conda:
        "envs/annotation.yaml"
     benchmark:
       "benchmarks/{sample}.annotateFiltered.benchmark.txt"
     shell:
        """
        export ANNOTSV={config[AnnotSVDir]}
        cpan YAML::XS
        cpan Sort::Key::Natural  
        {config[AnnotSVDir]}bin/AnnotSV -annotationsDir {config[AnnotSVDir]}share/AnnotSV/  -SvinputFile {input} -genomeBuild {config[genomeBuildGR]} outputDir {config[outPath]}results/ -outputFile {wildcards.sample}.annotatedFiltered 
        perl {config[knotAnnotSV]}knotAnnotSV.pl --configFile {config[knotAnnotSV]}config_insertions.yaml --annotSVfile {output[0]} --outDir {config[outPath]}results --genomeBuild hg19  
        """

rule annotateVerified:
     input:
       outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed"
     output:
       outPath + "results/{sample}.annotatedVerified.tsv",
       outPath + "results/{sample}.annotatedVerified.html"
     conda:
        "envs/annotation.yaml"
     benchmark:
       "benchmarks/{sample}.annotateVerified.benchmark.txt"
     shell:
        """
	export ANNOTSV={config[AnnotSVDir]}
        cpan YAML::XS
        cpan Sort::Key::Natural
        {config[AnnotSVDir]}bin/AnnotSV -annotationsDir {config[AnnotSVDir]}share/AnnotSV/  -SvinputFile {input} -genomeBuild {config[genomeBuildGR]} outputDir {config[outPath]}results/ -outputFile {wildcards.sample}.annotatedVerified 
        perl {config[knotAnnotSV]}knotAnnotSV.pl --configFile {config[knotAnnotSV]}config_insertions.yaml --annotSVfile {output[0]} --outDir {config[outPath]}results --genomeBuild hg19  
        """
