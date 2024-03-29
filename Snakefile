configfile: "config.yaml"

SAMPLES = ["BAM1_PREFIX","BAM2_PREFIX","BAM3_PREFIX"]
outPath = config["outPath"]
bamPath = config["bamPath"]
cramPath = config["cramPath"]


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

rule CramToBam:
    input:
        cram_file=cramPath + "{sample}.cram"
    output:
        bamPath + "{sample}.bam",
        bamPath + "{sample}.bam.bai"
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
        "retroseq.pl -discover -bam {input[0]} -output {output} -eref {config[HERVK_eref]} -id {params.identity}"
 
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
        "retroseq.pl -call -bam {input.bam} -input {input.discover} -ref {config[refHg19]} -output {output}"

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
    params:
        verificationLevel="low"
    conda:
       "envs/verification.yaml"
    log:
        "logs/call/{sample}.log"
    shell:
      """
      python {config[pythonScripts]}/assembleAndRepeatMasker.py {input[0]} {config[bamPath]}{wildcards.sample}.bam {config[outPath]} {config[RepeatMaskerPath]} {config[pythonScripts]} {config[element]} {params.verificationLevel} {output} 
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
       """ 
       bedtools window -w 500 -c -a {config[knownNR]} -b  {input}  > {output} 
       #rm {input}.temp.bed 
       """
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
	mkdir -p {config[outPath]}novel 
        bedtools sort -i  {input} >  {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        bedtools window -w 500 -v -a {config[outPath]}novel/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output}  
        rm {config[outPath]}novel/{wildcards.sample}.sorted.bed
        #rm {input}.temp.bed 
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
       """ 
       bedtools window -w 500 -c -a {config[knownNR]} -b  {input} > {output} 
       #rm {input}.temp.bed 
       """
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
	mkdir -p {config[outPath]}novel
        bedtools sort -i  {input} >  {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        bedtools window -w 500 -v -a {config[outPath]}novel/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output}  
        rm {config[outPath]}novel/{wildcards.sample}.sorted.bed 
        #rm {input}.temp.bed       
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
        sed "s/\\s*$/\\tMEI/" {input} > {input}.temp1.bed
        {config[AnnotSVDir]}bin/AnnotSV -annotationsDir {config[AnnotSVDir]}share/AnnotSV/  -SvinputFile {input}.temp1.bed -genomeBuild {config[genomeBuildGR]} -svtBEDcol 4 -outputDir {config[outPath]}results/ -outputFile {wildcards.sample}.annotatedFiltered 
        perl {config[knotAnnotSV]}knotAnnotSV.pl --configFile  envs/config_AnnotSV.yaml --annotSVfile {output[0]} --outDir {config[outPath]}results --genomeBuild hg19  
        rm {input}.temp1.bed
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
        sed "s/\\s*$/\\tMEI/" {input} > {input}.temp2.bed
        {config[AnnotSVDir]}bin/AnnotSV -annotationsDir {config[AnnotSVDir]}share/AnnotSV/  -SvinputFile {input}.temp2.bed -genomeBuild {config[genomeBuildGR]} -svtBEDcol 4 -outputDir {config[outPath]}results/ -outputFile {wildcards.sample}.annotatedVerified 
        perl {config[knotAnnotSV]}knotAnnotSV.pl --configFile  envs/config_AnnotSV.yaml --annotSVfile {output[0]} --outDir {config[outPath]}results --genomeBuild hg19  
        rm {input}.temp2.bed
	"""

