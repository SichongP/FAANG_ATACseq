
localrules: all, raw_multiqc, trim_multiqc, aggregateStats
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

#File name pattern
#{tissue}_{rep}_{assay}_{read}.fq.gz

rule all:
    input: 
        #expand("Results/mapping/stats/{tissue}_{rep}_{assay}.txt", tissue = TISSUES, rep = REPS, assay = ASSAYS),
        "Results/metrics/mapping_stats.csv",
        "Results/QC/raw/multiqc_report.html",
        "Results/QC/trimmed/multiqc_report.html",
#        expand("Results/figures/PEFragSize_{tissue}_{assay}.png", tissue = TISSUES, assay = ASSAYS),
#        "Results/figures/pearson_heatmap.png",
        expand("Results/deeptools/{tissue}_{rep}_{assay}.nonorm.bg", tissue = TISSUES, assay = ASSAYS, rep = REPS),
        "Results/figures/Fingerprint.png"

rule raw_qc:
    input: "raw_reads/{tissue}_{rep}_{assay}_{read}.fq.gz"
    output: "Results/QC/raw/{tissue}_{rep}_{assay}_{read}_fastqc.zip"
    conda: "env/fastqc.yaml"
    params: time=60
    threads: 1
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
    params: time=300
    threads: 1
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
    params: time=1200
    threads: 6
    shell:
     """
     bwa mem -t {threads} -R "@RG\\tID:{wildcards.tissue}_{wildcards.rep}_{wildcards.assay}\\tSM:{wildcards.tissue}_{wildcards.rep}_{wildcards.assay}\\tPL:illumina" -o {output} {config[ref]} {input}
     """

rule convert_sam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.sam"
    output: temp("Results/mapping/{tissue}_{rep}_{assay,\w+}.bam")
    conda: "env/samtools.yaml"
    params: time="120"
    threads: 4
    shell:
     """
     samtools view -bh -@ {threads} {input} > {output}
     """

rule sort_bam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.bam"
    output: temp("Results/mapping/{tissue}_{rep}_{assay}.sorted.bam")
    params: time="360"
    conda: "env/samtools.yaml"
    threads: 6
    shell:
     """
     cleanup() {{ rm -rf /scratch/pengsc/$SLURM_JOBID; }}
     trap cleanup EXIT
     mkdir -p /scratch/pengsc/$SLURM_JOBID/
     samtools sort -@ {threads} -m 1G -T /scratch/pengsc/$SLURM_JOBID/ -o {output} {input}
     """

rule index_bam:
    input: "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam"
    output: temp("Results/mapping/{tissue}_{rep}_{assay}.sorted.bam.bai")
    conda: "env/samtools.yaml"
    params: time="120"
    threads: 4
    shell:
     """
     samtools index -@ {threads} {input}
     """

rule markdup:
    input: bam = "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam", bai = "Results/mapping/{tissue}_{rep}_{assay}.sorted.bam.bai"
    output: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    conda: "env/sambamba.yaml"
    params: time="240"
    threads: 4
    shell:
     """
     MYTMPDIR=/scratch/pengsc/$SLURM_JOBID
     cleanup() {{ rm -rf $MYTMPDIR; }}
     trap cleanup EXIT
     mkdir -p $MYTMPDIR
     sambamba markdup --sort-buffer-size=8192 --io-buffer-size=1024 -t {threads} --tmpdir=$MYTMPDIR {input.bam} {output}
     """

rule bam_stat:
    input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}.txt"
    conda: "env/sambamba.yaml"
    params: time="50"
    threads: 4
    shell:
     """
     sambamba flagstat {input} > {output}
     """

rule aggregateStats:
    input: expand("Results/mapping/stats/{tissue}_{rep}_{assay}{ext}", tissue = TISSUES, rep = REPS, assay = ASSAYS, ext = [".txt", "_chrM.stats"])
    params: dir = "Results/mapping/stats/"
    output: "Results/metrics/mapping_stats.csv"
    conda: "env/stat_curator.yaml"
    script: "tools/flagstat_curator.py"

rule EstFragSize:
    input: expand("Results/mapping/{{tissue}}_{rep}_{{assay}}.markdup.sorted.bam", rep = REPS)
    output: fig = "Results/figures/PEFragSize_{tissue}_{assay}.png", stat = "Results/metrics/PEFragSize_{tissue}_{assay}.tsv", raw = "Results/metrics/FragSizeCounts_{tissue}_{assay}.txt"
    conda: "env/deeptools.yaml"
    params: time="240"
    threads: 8
    shell:
     """
     bamPEFragmentSize -hist {output.fig} -p {threads} -T "Fragment Size Distribution" -n 100000 --table {output.stat} --outRawFragmentLengths {output.raw} -b {input}
     """

rule bamCoverage_unfiltered:
    input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    output: "Results/deeptools/{tissue}_{rep}_{assay}.bigwig"
    conda: "env/deeptools.yaml"
    params: time="240"
    threads: 4
    shell:
     """
     bamCoverage -o {output} -of bigwig -p {threads} --effectiveGenomeSize {config[effective_size]} --normalizeUsing RPKM --exactScaling -b {input}
     """

rule makebedGraph_unfiltered:
    input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    output: "Results/deeptools/{tissue}_{rep}_{assay}.nonorm.bg"
    conda: "env/deeptools.yaml"
    params: time="240"
    threads: 4
    shell:
     """
     bamCoverage -o {output} -of bedgraph -p {threads} --effectiveGenomeSize {config[effective_size]} --exactScaling -b {input}
     """


rule multiBigWig_unfiltered:
    input: expand("Results/deeptools/{tissue}_{rep}_{assay}.bigwig", tissue = TISSUES, rep = REPS, assay = ASSAYS)
    output: "Results/deeptools/multiBigWig.npz"
    conda: "env/deeptools.yaml"
    threads: 8
    params: time="300"
    shell:
     """
     multiBigwigSummary bins -o {output} --smartLabels -bs 1000 -p {threads} -b {input}
     """
    
rule plotCorrelation_unfiltered:
    input: "Results/deeptools/multiBigWig.npz"
    output: fig = "Results/figures/pearson_heatmap.png", matrix = "Results/metrics/corrMatrix.tsv"
    conda: "env/deeptools.yaml"
    params: time = "240"
    threads: 6
    shell:
     """
     plotCorrelation --corData {input} -c pearson -p heatmap -o {output.fig} --skipZeros --removeOutliers --outFileCorMatrix {output.matrix}
     """

rule getMitoReads:
    input: "Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam"
    output: "Results/mapping/{tissue}_{rep}_{assay}.chrM.bam"
    conda: "env/samtools.yaml"
    params: time = "240"
    threads: 4
    shell:
     """
     samtools view -bh -@ {threads} {input} chrM > {output}
     """    

rule chrMStats:
    input: "Results/mapping/{tissue}_{rep}_{assay}.chrM.bam"
    output: "Results/mapping/stats/{tissue}_{rep}_{assay}_chrM.stats"
    conda: "env/samtools.yaml"
    params: time = "240"
    threads: 4
    shell:
     """
     samtools flagstat -@ {threads} {input} > {output}
     """

rule plotFingerprint_all:
    input:
        bam = expand("Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam", rep = REPS, tissue = TISSUES, assay = ASSAYS),
        bai = expand("Results/mapping/{tissue}_{rep}_{assay}.markdup.sorted.bam.bai", rep = REPS, tissue = TISSUES, assay = ASSAYS)
    output:
        plot = "Results/figures/Fingerprint.png", metrics = "Results/metrics/Fingerprint.txt"
    conda: "env/deeptools.yaml"
    params: time = "360"
    threads: 8
    log: "Results/logs/deeptools/fingerprint.log"
    shell:
     """
     plotFingerprint --bamfiles {input.bam} -o {output.plot} -p {threads} --skipZeros -T "Fingerprint plot" --smartLabels --outQualityMetrics {output.metrics}
     """

