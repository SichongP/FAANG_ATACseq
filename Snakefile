
localrules: all, raw_multiqc, trim_multiqc, aggregateStats
ruleorder: markdup > sort_bam > subSample > mapping
import re, os

configfile: "config.yaml"

FILES = []
for dirpath, dirname, filenames, in os.walk("raw_reads/"):
    FILES.extend(filenames)
    break

BASENAMES = []
for filename in FILES:
    name = re.match(r"(.+?)\.fq\.gz", filename)
    if name:
        BASENAMES.append(name.group(1))

TISSUES = config['tissues']
REPS = config['reps']
ASSAYS  = config['assays']
READS = config['reads']
REF = config['ref']
EFFECTIVE_SIZE = config['effective_size']


def get_fraction(wildcards):
    assay = wildcards.assay
    tissue = wildcards.tissue
    rep = wildcards.rep
    target = wildcards.count
    with open("Results/mapping/stats/{}_{}_{}.txt".format(tissue, rep, assay)) as stat_file:
        line = stat_file.readlines()[0]
#        print(line)
#        print(line.split('+')[0].strip())
        total = float(line.split('+')[0].strip())
    return min(round(int(target) * 1000000 / total, 2), 1)


#File name pattern
#{tissue}_{rep}_{assay}_{read}.fq.gz

rule all:
    input: 
        #expand("Results/mapping/stats/{tissue}_{rep}_{assay}.txt", tissue = TISSUES, rep = REPS, assay = ASSAYS),
        "Results/metrics/mapping_stats.csv",
        "Results/metrics/mapping_stats_Liver.csv",
        "Results/QC/raw/multiqc_report.html",
        "Results/QC/trimmed/multiqc_report.html",
        expand("Results/figures/PEFragSize_{tissue}_{assay}{count}.png", tissue = TISSUES, assay = ASSAYS, count = ['', '_40M', '_50M', '_60M']),
        "Results/figures/pearson_heatmap.png",
        expand("Results/deeptools/{tissue}_{rep}_{assay}{count}.nonorm.bg", tissue = TISSUES, assay = ASSAYS, rep = REPS, count = ['', '_40M', '_50M', '_60M']),
        "Results/figures/Fingerprint.png",
        expand("Results/figures/Fingerprint_{tissue}.png", tissue = TISSUES),
#        expand("Results/mapping/{tissue}_{rep}_{assay}.shifted.sorted.bam", tissue = TISSUES, assay = ASSAYS, rep = REPS),
#        expand("Results/mapping/{tissue}_{rep}_{assay}.shifted.sorted.bam.bai", tissue = TISSUES, assay = ASSAYS, rep = REPS),
        expand("Results/mapping/{tissue}_{rep}_{assay}.dedup.bam", tissue = TISSUES, assay = ASSAYS, rep = REPS),
        expand("Results/mapping/Liver_{rep}_NucleiFirst.dedup.bam", rep = REPS),
        expand("Results/mapping/Liver_{rep}_NucleiFirst.dedup.bam.bai", rep = REPS),
        expand("Results/mapping/{tissue}_{rep}_{assay}.dedup.bam.bai", tissue = TISSUES, assay = ASSAYS, rep = REPS),
        expand("Results/mapping/{tissue}_{rep}_{assay}.chrM.bam", tissue = TISSUES, assay = ASSAYS, rep = REPS),
        expand("Results/mapping/Liver_{rep}_NucleiFirst.chrM.bam", rep = REPS),
        expand("Results/hmmratac/{tissue}_{rep}_{assay}{count}_peaks.gappedPeak", tissue = TISSUES, assay = ASSAYS, rep = REPS, count = ['', '_40M', '_50M', '_60M']),
#        expand("Results/hmmratac/Liver_{rep}_NucleiFirst_peaks.gappedPeak", rep = REPS),
#        expand("Results/macs2/{tissue}_{rep}_{assay}{count}_peaks.broadPeak", tissue = TISSUES, assay = ASSAYS, rep = REPS, count = ['', '_40M', '_50M', '_60M']),
#        expand("Results/macs2/Liver_{rep}_NucleiFirst_peaks.broadPeak", rep = REPS),
        expand("Results/figures/{tissue}_{rep}_{assay}{count}.enrichment{program}.png", tissue = TISSUES, assay = ASSAYS, rep = REPS, program = ['','.macs2'], count = ['', '_40M', '_50M', '_60M']),
        expand("Results/figures/Liver_{rep}_NucleiFirst.enrichment{program}.png", rep = REPS, program = ['','.macs2']),
        "Results/figures/PCA.png",
        expand("Results/figures/PCA_{tissue}.png", tissue = TISSUES),
        expand("Results/figures/{tissue}_{rep}_{assay}_Fingerprint.png", tissue = TISSUES, assay = ASSAYS, rep = REPS)


rule raw_qc:
    input: "raw_reads/{tissue}_{rep}_{assay}_{read}.fq.gz"
    output: "Results/QC/raw/{tissue}_{rep}_{assay}_{read}_fastqc.zip"
    conda: "env/fastqc.yaml"
    resources: mem_mb = 500, time_min = 60
    params: partition = 'low2'
    shell:
     """
     fastqc -o Results/QC/raw/ {input}
     """

rule raw_multiqc:
    input: expand("Results/QC/raw/{tissue}_{rep}_{assay}_{read}_fastqc.zip", tissue = TISSUES, rep = REPS, assay = ASSAYS, read = READS)
    output: "Results/QC/raw/multiqc_report.html"
    conda: "env/fastqc.yaml"
    shell:
     """
     multiqc -o Results/QC/raw/ Results/QC/raw/
     """

rule trim:
    input: expand("raw_reads/{{tissue}}_{{rep}}_{{assay}}_{read}.fq.gz", read = READS)
    output: temp(expand("Results/trimming/{{tissue}}_{{rep}}_{{assay}}_{read}.fq.gz", read = ["val_1", "val_2"])), expand("Results/QC/trimmed/{{tissue}}_{{rep}}_{{assay}}_{read}_fastqc.zip", read = ["val_1", "val_2"])
    conda: "env/trim_galore.yaml"
    params: partition = 'med2'
    resources: time_min=1000, mem_mb=8000, mem_mb_med=8000, cpus_med=1
    shell:
     """
     trim_galore --paired -o Results/trimming/ --fastqc_args "-o Results/QC/trimmed/" --basename {wildcards.tissue}_{wildcards.rep}_{wildcards.assay} {input}
     """

rule trim_multiqc:
    input: expand("Results/QC/trimmed/{tissue}_{rep}_{assay}_{read}_fastqc.zip", tissue = TISSUES, rep = REPS, assay = ASSAYS, read = ["val_1", "val_2"])
    output: "Results/QC/trimmed/multiqc_report.html"
    conda: "env/fastqc.yaml"
    shell:
     """
     multiqc -o Results/QC/trimmed/ Results/QC/trimmed/
     """

rule mapping:
    input: expand("Results/trimming/{{tissue}}_{{rep}}_{{assay}}_{read}.fq.gz", read = ["val_1", "val_2"])
    output: temp("Results/mapping/{tissue}_{rep}_{assay,\w+}.sam")
    conda: "env/bwa.yaml"
    params: partition = 'bmm'
    resources: time_min=5000, cpus=10, cpus_bmm=10, mem_mb=lambda wildcards, attempt: 15000+5000*attempt, mem_mb_bmm=lambda wildcards, attempt: 15000+5000*attempt
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.tissue}_{wildcards.rep}_{wildcards.assay}\\tSM:{wildcards.tissue}_{wildcards.rep}_{wildcards.assay}\\tPL:illumina" -o {output} {config[ref]} {input}
     """

rule convert_sam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.sam"
    output: temp("Results/mapping/{tissue}_{rep}_{assay,\w+}.bam")
    conda: "env/samtools.yaml"
    params: partition = 'med2'
    resources: cpus=4, mem_mb=4000, cpus_med=4, mem_mb_med=4000, time_min = 500
    shell:
     """
     samtools view -bh -@ {resources.cpus} {input} > {output}
     """

rule sort_bam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.bam"
    output: temp("Results/mapping/{tissue}_{rep}_{assay}.sorted.bam")
    params: partition = "med2"
    conda: "env/samtools.yaml"
    resources: time_min=600, cpus=6, cpus_med = 6, mem_mb=lambda wildcards, attempt: 6000+2000*attempt, mem_mb_med=lambda wildcards, attempt: 6000+2000*attempt
    shell:
     """
     cleanup() {{ rm -rf /scratch/pengsc/$SLURM_JOBID; }}
     trap cleanup EXIT
     mkdir -p /scratch/pengsc/$SLURM_JOBID/
     samtools sort -@ {resources.cpus} -m 1G -T /scratch/pengsc/$SLURM_JOBID/ -o {output} {input}
     """

rule index_bam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam"
    output: "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam.bai"
    conda: "env/samtools.yaml"
    params: partition = 'low2'
    resources: cpus=4, mem_mb=4000, time_min = 600
    shell:
     """
     samtools index -@ {resources.cpus} {input}
     """

rule markdup:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam.bai"
    output: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    conda: "env/sambamba.yaml"
    params: partition = 'bmm'
    resources: time_min=600, cpus=6, cpus_bmm=6, mem_mb=lambda wildcards, attempt: 15000+10000*attempt, mem_mb_bmm=lambda wildcards, attempt: 15000+10000*attempt
    shell:
     """
     MYTMPDIR=/scratch/pengsc/$SLURM_JOBID
     cleanup() {{ rm -rf $MYTMPDIR; }}
     trap cleanup EXIT
     mkdir -p $MYTMPDIR
     sambamba markdup --sort-buffer-size=8192 --io-buffer-size=1024 -t {resources.cpus} --tmpdir=$MYTMPDIR {input.bam} {output}
     """

rule subSample:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam.bai", stat = "Results/mapping/stats/{tissue}_{rep}_{assay}.txt"
    output: bam = "Results/mapping/{tissue}_{rep}_{assay}_{count}M.markdup.sorted.bam"
#, bai = "Results/mapping/{tissue}_{rep}_{assay}_{count}M.markdup.sorted.bam.bai"
    conda: "env/sambamba.yaml"
    params: partition = 'med2', frac = get_fraction
    resources: cpus=4, cpus_med=4, time_min=900, mem_mb = 10000, mem_mb_med = 10000
    shell:
     """
     sambamba view -s {params.frac} -t {resources.cpus} -o {output.bam} -f bam {input.bam}
     """

#rule index_markdup:
 #   input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
  #  output: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam.bai"
   # conda: "env/samtools.yaml"
   # params: time="120"
   # threads: 4
   # shell:
   #  """ 
   #  samtools index -@ {threads} {input}
   #  """

rule bam_stat:
    input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}.txt"
    conda: "env/sambamba.yaml"
    params: partition = 'low2'
    resources: cpus=4, mem_mb=8000
    shell:
     """
     sambamba flagstat {input} > {output}
     """

rule aggregateStats:
    input: expand("Results/mapping/stats/{tissue}_{rep}_{assay}{count}{ext}", tissue = TISSUES, rep = REPS, assay = ASSAYS, ext = [".txt", "_chrM.stats", "_filtered.txt", "_dedup.txt"], count = ['', '_40M', '_50M', '_60M']), expand("Results/mapping/stats/Liver_{rep}_NucleiFirst{ext}", rep = REPS, ext = [".txt", "_chrM.stats", "_filtered.txt", "_dedup.txt"])
    params: dir = "Results/mapping/stats/"
    output: "Results/metrics/mapping_stats.csv"
    conda: "env/stat_curator.yaml"
    script: "tools/flagstat_curator.py"

rule aggregateStats_Liver:
    input: expand("Results/mapping/stats/Liver_{rep}_{assay}{count}{ext}", rep = REPS, assay = ['Nuclei', 'Tissue'], ext = [".txt", "_chrM.stats", "_filtered.txt", "_dedup.txt"], count = ['', '_5M', '_10M', '_15M', '_20M', '_25M', '_30M', '_35M', '_40M', '_50M', '_55M', '_60M', '_65M', '_70M', '_75M', '_80M', '_85M', '_90M', '_95M', '_100M', '_105M', '_110M', '_115M', '_120M', '_125M', '_130M']), expand("Results/mapping/stats/Liver_{rep}_NucleiFirst{ext}", rep = REPS, ext = [".txt", "_chrM.stats", "_filtered.txt", "_dedup.txt"])
    params: dir = "Results/mapping/stats/"
    output: "Results/metrics/mapping_stats_Liver.csv"
    conda: "env/stat_curator.yaml"
    script: "tools/flagstat_curator.py"


rule EstFragSize:
    input: expand("Results/mapping/{{tissue}}_{rep}_{{assay}}.markdup.sorted.bam", rep = REPS)
    output: fig = "Results/figures/PEFragSize_{tissue}_{assay}.png", stat = "Results/metrics/PEFragSize_{tissue}_{assay}.tsv", raw = "Results/metrics/FragSizeCounts_{tissue}_{assay}.txt"
    wildcard_constraints: tissue="[a-zA-Z]+"
    conda: "env/deeptools.yaml"
    params: partition='low2'
    resources: cpus=8, mem_mb = 4000, time_min = 300
    shell:
     """
     bamPEFragmentSize -hist {output.fig} -p {resources.cpus} -T "Fragment Size Distribution" -n 100000 --table {output.stat} --outRawFragmentLengths {output.raw} -b {input}
     """

rule bamCoverage:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam"
    output: "Results/deeptools/{tissue}_{rep}_{assay}.bigwig"
    conda: "env/deeptools.yaml"
    params: partition = 'low2'
    resources: time_min = 640, cpus = 4, mem_mb = 8000
    threads: 4
    shell:
     """
     bamCoverage -o {output} -of bigwig -p {resources.cpus} --effectiveGenomeSize {config[effective_size]} --normalizeUsing RPKM --exactScaling -b {input.bam}
     """

rule makebedGraph:
    input: "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam"
    output: "Results/deeptools/{tissue}_{rep}_{assay}.nonorm.bg"
    conda: "env/deeptools.yaml"
    params: partition = 'low2'
    resources: cpus=4, mem_mb=8000, time_min=300
    shell:
     """
     bamCoverage -o {output} -of bedgraph -p {resources.cpus} --effectiveGenomeSize {config[effective_size]} --exactScaling -b {input}
     """


rule multiBigWig:
    input: expand("Results/deeptools/{tissue}_{rep}_{assay}{count}.bigwig", tissue = TISSUES, rep = REPS, assay = ASSAYS, count=['', '_40M', '_50M', '_60M'])
    output: "Results/deeptools/multiBigWig.npz"
    conda: "env/deeptools.yaml"
    resources: cpus=8, time_min=900, mem_mb=20000, mem_mb_med=20000, cpus_med=8
    params: partition='med2'
    shell:
     """
     multiBigwigSummary bins -o {output} --smartLabels -bs 1000 -p {resources.cpus} -b {input}
     """

rule multiBigWig_byTissue:
    input: expand("Results/deeptools/{{tissue}}_{rep}_{assay}{count}.bigwig", rep = REPS, assay = ASSAYS, count=['', '_40M', '_50M', '_60M'])
    output: "Results/deeptools/{tissue}_multiBigWig.npz"
    conda: "env/deeptools.yaml"
    resources: cpus=8, time_min=900, mem_mb=20000, mem_mb_med=20000, cpus_med=8
    params: partition='med2'
    shell:
     """
     multiBigwigSummary bins -o {output} --smartLabels -bs 1000 -p {resources.cpus} -b {input}
     """    

rule plotCorrelation_unfiltered:
    input: "Results/deeptools/multiBigWig.npz"
    output: fig = "Results/figures/pearson_heatmap.png", matrix = "Results/metrics/corrMatrix.tsv"
    conda: "env/deeptools.yaml"
    params: partition='low2'
    resources: cpus=1, time_min=240, mem_mb=8000
    shell:
     """
     plotCorrelation --corData {input} -c pearson -p heatmap -o {output.fig} --skipZeros --removeOutliers --outFileCorMatrix {output.matrix}
     """

rule getMitoReads:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam.bai"
    output: "Results/mapping/{tissue}_{rep}_{assay}.chrM.bam"
    conda: "env/samtools.yaml"
    params: partition='low2'
    resources: cpus=4, time_min=240, mem_mb=10000
    shell:
     """
     samtools view -bh -@ {resources.cpus} {input.bam} chrM > {output}
     """    

rule chrMStats:
    input: "Results/mapping/{tissue}_{rep}_{assay}.chrM.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}_chrM.stats"
    conda: "env/samtools.yaml"
    params: partition='low2'
    resources: cpus=4, time_min=240, mem_mb=4000
    shell:
     """
     samtools flagstat -@ {resources.cpus} {input} > {output}
     """

rule plotFingerprint_all:
    input:
        bam = expand("Results/mapping/{tissue}_{rep}_{assay}.dedup.bam", rep = REPS, tissue = TISSUES, assay = ASSAYS)
    output:
        plot = "Results/figures/Fingerprint.png", metrics = "Results/metrics/Fingerprint.txt"
    conda: "env/deeptools.yaml"
    params: partition='bmm'
    resources: cpus=8, cpus_bmm=8, mem_mb=20000, mem_mb_bmm=20000, time_min=360
#    log: "Results/logs/deeptools/fingerprint.log"
    shell:
     """
     plotFingerprint --bamfiles {input.bam} -o {output.plot} -p {resources.cpus} --skipZeros -T "Fingerprint plot" --smartLabels --outQualityMetrics {output.metrics} --extendReads
     """

rule plotFingerprint_Tissue:
    input:
        bam = expand("Results/mapping/{{tissue}}_{rep}_{assay}.dedup.bam", rep = REPS, assay = ASSAYS)
    output:
        plot = "Results/figures/Fingerprint_{tissue}.png", metrics = "Results/metrics/Fingerprint_{tissue}.txt"
    conda: "env/deeptools.yaml"
    params: partition='bmm'
    resources: cpus=8, cpus_bmm=8, mem_mb=20000, mem_mb_bmm=20000, time_min=360
    shell:
     """
     plotFingerprint --bamfiles {input.bam} -o {output.plot} -p {resources.cpus} --skipZeros -T "Fingerprint plot" --smartLabels --outQualityMetrics {output.metrics} --extendReads
     """

rule plotFingerprint_Lib:
    input:
        bam = expand("Results/mapping/{{tissue}}_{{rep}}_{{assay}}{count}.dedup.bam", count = ['', '_40M', '_50M', '_60M'])
    output:
        plot = "Results/figures/{tissue}_{rep}_{assay}_Fingerprint.png", metrics = "Results/metrics/{tissue}_{rep}_{assay}_Fingerprint.txt"
    conda: "env/deeptools.yaml"
    params: partition='bmm'
    resources: cpus=8, cpus_bmm=8, mem_mb=20000, mem_mb_bmm=20000, time_min=360
    shell:
     """
     plotFingerprint --bamfiles {input.bam} -o {output.plot} -p {resources.cpus} --skipZeros -T "Fingerprint plot" --smartLabels --outQualityMetrics {output.metrics} --extendReads
     """


rule filterBam:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam.bai"
    output: bam = "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam.bai"
    params: partition='low2'
    resources: cpus=4, time_min=360, mem_mb=10000
    conda: "env/samtools.yaml"
    shell:
     """
     samtools idxstats {input.bam} | cut -f 1 | grep -v "chrUn" | grep -v "chrM" | xargs samtools view -bh -@ {resources.cpus} -F 2820 {input.bam} > {output.bam}
     samtools index -@ {resources.cpus} {output.bam}
     """

rule fileterDup:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam.bai"
    output: bam = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam.bai"
    params: partition='low2'
    resources: cpus=4, time_min=360, mem_mb=8000
    conda: "env/samtools.yaml"
    shell:
     """
     xargs samtools view -bh -@ {resources.cpus} -F 1024 {input.bam} > {output.bam}
     samtools index -@ {resources.cpus} {output.bam}
     """

rule filtered_stat:
    input: "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}_filtered.txt"
    conda: "env/sambamba.yaml"
    params: partition='low2'
    resources: cpus=4, time_min=50, mem_mb=4000
    shell:
     """
     sambamba flagstat {input} > {output}
     """

rule dedupStat:
    input: "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}_dedup.txt"
    conda: "env/sambamba.yaml"
    params: partition='low2'
    resources: cpus=4, time_min=50, mem_mb=4000
    shell:
     """
     sambamba flagstat {input} > {output}
     """


rule shiftReads:
    input: "Results/mapping/{tissue}_{rep}_{assay}.filtered.bam"
    output: "Results/mapping/{tissue}_{rep}_{assay}.shifted.bam"
    conda: "env/deeptools.yaml"
    params: partition = 'low2'
    resources: cpus=4, time_min=360, mem_mb=8000
    shell:
     """
     alignmentSieve -b {input} -o {output} -p {resources.cpus} --ATACshift 
     """

rule hmmratac:
    input: 
        bam = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam",
        bai = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam.bai",
        geno_info = "genome.info"
    output: "Results/hmmratac/{tissue}_{rep}_{assay}_peaks.gappedPeak"
    conda: "env/hmmratac.yaml"
    params: partition='bmm'
    resources: cpus=1, cpus_bmm=1, mem_mb=80000, mem_mb_bmm=80000, time_min=1000
    shell:
     """
     java -jar /home/pengsc/bin/miniconda3/envs/hmmratac/share/hmmratac-1.2.9-0/HMMRATAC.jar -Xms30g -b {input.bam} --peaks True -i {input.bai} --bedgraph True -g {input.geno_info} -o Results/hmmratac/{wildcards.tissue}_{wildcards.rep}_{wildcards.assay} --threshold 2 --score fc
     """

rule macs2:
    input:
        bam = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam",
        bai = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam.bai"
    output: "Results/macs2/{tissue}_{rep}_{assay}_peaks.broadPeak", "Results/macs2/{tissue}_{rep}_{assay}_treat_pileup.bdg", "Results/macs2/{tissue}_{rep}_{assay}_control_lambda.bdg", "Results/macs2/{tissue}_{rep}_{assay}_peaks.gappedPeak"
    conda: "env/macs2.yaml"
    params: partition='bmm'
    resources: cpus=1, time_min=600, cpus_bmm=1, mem_mb=10000, mem_mb_bmm=10000
    shell:
     """
     macs2 callpeak -t {input.bam} -f BAMPE --keep-dup all --outdir Results/macs2 -n {wildcards.tissue}_{wildcards.rep}_{wildcards.assay} -B --SPMR --cutoff-analysis --broad --nolambda
     """

rule plotPCA:
    input: "Results/deeptools/multiBigWig.npz"
    output: "Results/figures/PCA.png"
    conda: "env/deeptools.yaml"
    params: partition='bml'
    resources: cpus=1, time_min=300, mem_mb=6000
    shell:
     """
     plotPCA -in {input} -o {output} --outFileNameData Results/figures/PCA_log2.tab --transpose --log2 -T "PCA plot of read coverage" 
     """

rule plotPCA_byTissue:
    input: "Results/deeptools/{tissue}_multiBigWig.npz"
    output: "Results/figures/PCA_{tissue}.png"
    conda: "env/deeptools.yaml"
    params: partition='bml'
    resources: cpus=1, time_min=300, mem_mb=6000
    shell:
     """ 
     plotPCA -in {input} -o {output} --outFileNameData Results/figures/PCA_log2_{wildcards.tissue}.tab --transpose --log2 -T "PCA plot of read coverage" 
     """

rule plotEnrichment:
    input: bams = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam", bed = "Results/hmmratac/{tissue}_{rep}_{assay}_peaks.gappedPeak"
    output: png = "Results/figures/{tissue}_{rep}_{assay}.enrichment.png", tab = "Results/figures/{tissue}_{rep}_{assay}.enrichment.tab"
    conda: "env/deeptools.yaml"
    params: partition='med2'
    resources: cpus=4, cpus_med=4, mem_mb=10000, mem_mb_med=10000, time_min=2400
    shell:
     """
     plotEnrichment -b {input.bams} --BED {input.bed} -o {output.png} --smartLabels -T "Fraction of Reads in Peaks (FRiP), {wildcards.tissue} {wildcards.rep} {wildcards.assay}" --outRawCounts {output.tab} -p {resources.cpus}
     """

rule plotEnrichment_macs2:
    input: bams = "Results/mapping/{tissue}_{rep}_{assay}.dedup.bam", bed = "Results/macs2/{tissue}_{rep}_{assay}_peaks.gappedPeak"
    output: png = "Results/figures/{tissue}_{rep}_{assay}.enrichment.macs2.png", tab = "Results/figures/{tissue}_{rep}_{assay}.enrichment.macs2.tab"
    conda: "env/deeptools.yaml"
    params: partition='med2'
    resources: cpus=4, cpus_med=4, mem_mb=10000, mem_mb_med=10000, time_min=2400
    shell:
     """
     plotEnrichment -b {input.bams} --BED {input.bed} -o {output.png} --smartLabels -T "Fraction of Reads in Peaks (FRiP), {wildcards.tissue} {wildcards.rep} {wildcards.assay}" --outRawCounts {output.tab} -p {threads}
     """

