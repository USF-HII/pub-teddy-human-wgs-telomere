# Snakefile

shell.prefix("set -euo pipefail; set -x;")

from os.path import abspath, basename, dirname, exists, join, splitext

BASE, RUN, WORK = os.environ["BASE"], os.environ["RUN"], os.environ["WORK"]

SINGULARITY_CMD = " ".join(["/shares/hii/sw/singularity/2.5.1/bin/singularity", "exec", "--bind /shares:/shares", "--bind /hii/work:/hii/work"])

SINGULARITY_IMG = os.environ.get("HII_SINGULARITY_IMAGE")

SINGULARITY = " ".join([SINGULARITY_CMD, SINGULARITY_IMG])

RSCRIPT = " ".join([SINGULARITY_CMD, os.environ.get("RSCRIPT_IMAGE")])

INPUT_FASTQ_DIR = os.environ.get("INPUT_FASTQ_DIR")

INPUT_DIR_BAM = {
    "grch37": join(os.environ["OUTPUT_TEDDY"], "grch37/bams"),
    "grch38": join(os.environ["OUTPUT_TEDDY"], "grch38/bams")}

REF_GENOME = {
    "grch37": "ref/gotcloud/hs37d5-db142-v1/hs37d5.fa",
    "grch38": "ref/gotcloud/hs38DH-db142-v1/hs38DH.fa"}

if os.environ.get("COMBINE_METHODS"):
    COMBINE_METHODS = os.environ["COMBINE_METHODS"].split(",")
else:
    COMBINE_METHODS = [
        "computel_merge.tsv",
        "motif_counter.tsv",
        "qmotif_merge.tsv",
        "telomerecat_telbam2length.csv",
        "telseq_merge.tsv"
        ]

WGS_METRICS_GRCH37 = join(os.environ.get("OUTPUT_TEDDY"),
    "grch37/picard_metrics_merge/wgs_metrics.tsv")

WGS_METRICS_GRCH38 = join(os.environ.get("OUTPUT_TEDDY"),
    "grch38/picard_metrics_merge/wgs_metrics.tsv")

if os.environ.get("SAMPLES"):
    SAMPLES = os.environ.get("SAMPLES").split(",")
else:
    SAMPLES, = glob_wildcards(join(INPUT_DIR_BAM["grch38"], "{sample}.bam"))

if RUN == "test":
    SAMPLES = SAMPLES[:2]

elif RUN == "test-100":
    SAMPLES = SAMPLES[:100]

#=====================================================================================================
# Targets
#=====================================================================================================

rule all:
    input:
        join(WORK, "computel_merge.tsv"),
        join(WORK, "motif_counter.tsv"),
        join(WORK, "qmotif_merge.tsv"),
        join(WORK, "telomerecat_telbam2length.csv"),
        join(WORK, "telseq_merge.tsv"),
        expand(join(WORK, "correlation/{fname}"),
               fname=["pearson.pdf", "spearman.pdf", "pearson.txt", "spearman.txt", "telomere.txt"])

rule target_combine_methods:
    input:
        join(WORK, "combined_methods.tsv")

rule target_correlation:
    input:
        expand(join(WORK, "correlation/{fname}"),
               fname=["pearson.pdf", "spearman.pdf", "pearson.txt", "spearman.txt", "telomere.txt"])

rule target_telomerecat:
    input:
        join(WORK, "telomerecat_telbam2length.csv")

rule target_telseq_merge:
    input:
        join(WORK, "telseq_merge.tsv")

rule target_computel:
    input:
        join(WORK, "computel_merge.tsv")

rule target_qmotif:
    input:
        join(WORK, "qmotif_merge.tsv")

rule target_motif_counter:
    input:
        join(WORK, "motif_counter.tsv")

#=====================================================================================================
# Functions
#=====================================================================================================
def get_coverage(fname):
    """
    Reads a single merged picard_metrics_merge/wgs_metrics.tsv file (from $OUTPUT_TEDDY)
    and returns a dictionary with (sample, mean_coverage) key pairs.
    """

    coverage = {}

    with open(fname) as f:
        header = next(f).rstrip("\n").split("\t")
        sample_column_idx = header.index("SAMPLE")
        mean_coverage_idx = header.index("MEAN_COVERAGE")

        for line in f:
            line = line.rstrip("\n").split("\t")
            coverage[line[sample_column_idx]] = float(line[mean_coverage_idx])

    return coverage

#=====================================================================================================
# Rules
#=====================================================================================================

localrules: combine_methods, telseq_merge, computel_merge, qmotif_merge, motif_counter_calc

#=====================================================================================================
# Telomerecat
#=====================================================================================================

rule telomerecat:
    input:
        bam=join(INPUT_DIR_BAM["grch38"], "{sample}.bam")
    output:
        telbam=join(WORK, "telomerecat/{sample}_telbam.bam")
    params:
        job_name="{sample}-telomerecat", log="telomerecat/{sample}.log", cpus="1", mem="6G", time="8-0",
        output_prefix=join(WORK, "telomerecat/{sample}")
    shadow:
        "shallow"
    shell:
        """
        {SINGULARITY} telomerecat bam2telbam \
            -v 2 \
            {input.bam}

        mv {wildcards.sample}_telbam.bam {output.telbam}
        """

rule telomerecat_telbam2length:
    input:
        telbams=expand(join(WORK, "telomerecat/{sample}_telbam.bam"), sample=SAMPLES)
    output:
        csv=join(WORK, "telomerecat_telbam2length.csv")
    params:
        job_name="telomerecat_telbam2length", log="telomerecat_telbam2length.log", cpus="8", mem="60G", time="2-0"
    shadow:
        "shallow"
    run:
        for telbam in input.telbams:
            shell(f"/bin/ln -s {telbam}")

        shell("{SINGULARITY} telomerecat telbam2length -v 2 -e --output {output.csv} *.bam")
        shell("{SINGULARITY} /bin/sed -i 's/\.bam,/,/' {output.csv}")

#=====================================================================================================
# Telseq
#=====================================================================================================

rule telseq:
    input:
        bam=join(INPUT_DIR_BAM["grch38"], "{sample}.bam")
    output:
        tsv=join(WORK, "telseq/{sample}.tsv")
    params:
        job_name="{sample}-telseq", log="telseq/{sample}.log", cpus="2", mem="12G", time="1-0"
    shell:
        """
        {SINGULARITY} telseq -u -r 150 -k 12 {input.bam} > {output.tsv}
        """

rule telseq_merge:
    input:
        tsvs=expand(join(WORK, "telseq/{sample}.tsv"), sample=SAMPLES)
    output:
        tsv=join(WORK, "telseq_merge.tsv")
    params:
        job_name="telseq_merge", log="telseq_merge.log", cpus="1", mem="10G", time="0-1"
    run:
        with open(input.tsvs[0]) as f_in:
            header = f_in.readlines()[2].rstrip("\n").split("\t")

        with open(output.tsv, "w") as f_out:
            print("\t".join(header[2:]), file=f_out)

            for fname in input.tsvs:
                with open(fname) as f_in:
                    sample = splitext(basename(fname))[0]
                    line = f_in.readlines()[3].rstrip("\n").split("\t")
                    print("\t".join([sample] + line[3:]), file=f_out)

#=====================================================================================================
# Computel
#=====================================================================================================

rule computel:
    input:
        fastq_1=join(INPUT_FASTQ_DIR, "{sample}_R1.fastq.gz"),
        fastq_2=join(INPUT_FASTQ_DIR, "{sample}_R2.fastq.gz")
    output:
        join(WORK, "computel/{sample}/tel.length.xls")
    params:
        job_name="computel-{sample}", log="computel/{sample}.log", cpus="4", mem="10G", time="1-0",
        output_prefix=join(WORK, "computel/{sample}")
    shell:
        """
        mkdir -p {params.output_prefix}

        {SINGULARITY} bash /computel/computel.sh \
            -1 {input.fastq_1} \
            -2 {input.fastq_2} \
            -pattern TTAGGG \
            -nchr 23 \
            -minseed 54 \
            -o {params.output_prefix}
        """

rule computel_merge:
    """
    Parse the telomere length value from each tel.length.xls file into a single file
    with the format: <sample>\t<telomere_length>
    """

    """
    Here is the example format of each tel.length.xls file (tab-separated) with line numbers (n|) prefixed:

        0|coverage.file   <path...>/MICH245849050730_P1/align/tel.align.bam.coverage.txt
        1|reads           <path...>/MICH245849050730_P1_R1.fastq.gz, <path...>/MICH245849050730_P1_R2.fastq.gz
        2|read.length     151
        3|pattern.length  c(pattern = 6)
        4|base.cov        42.3044334801409
        5|num.haploid.chr 23
        6|tel.length      c(pattern = 5226.75971537241)
        7|enome.length    3244610000
        8|min.seed        54
    """

    input:
        tsvs=expand(join(WORK, "computel/{sample}/tel.length.xls"), sample=SAMPLES)
    output:
        tsv=join(WORK, "computel_merge.tsv")
    params:
        job_name="computel_merge", log="computel_merge.log", cpus="1", mem="10G", time="0-1"
    run:
        header = ["Sample", "Tel_Length"]

        with open(output.tsv, "w") as f_out:
            print("\t".join(header), file=f_out)

            for fname in input.tsvs:
                sample = dirname(fname).split("/")[-1]

                with open(fname) as f_in:
                    length = f_in.readlines()[6].rstrip("\n").split(" = ")[1][:-1]
                    print("\t".join([sample, length]), file=f_out)


#=====================================================================================================
# qMotif
#=====================================================================================================

rule qmotif:
    """
    Tool only works for GRCh37 bam files
    """
    input:
        bam=join(INPUT_DIR_BAM["grch37"], "{sample}.bam")
    output:
        join(WORK, "qmotif/{sample}/{sample}.xml")
    params:
        job_name="qmotif-{sample}", log="qmotif/{sample}.log", cpus="1", mem="10G", time="0-1",
        output_prefix=join(WORK, "qmotif/{sample}")
    shell:
        """
        mkdir -p {params.output_prefix}

        {SINGULARITY} java -Xmx20g -jar /usr/local/bin/qmotif-1.2.jar \
            -n 4 \
            --bam {input.bam} \
            --bai {input.bam}.bai \
            --log {params.output_prefix}/out.log \
            --loglevel INFO \
            -ini {BASE}/metadata/qmotif.ini \
            -o {params.output_prefix}/{wildcards.sample}.xml \
            -o {params.output_prefix}/{wildcards.sample}.bam
        """

rule qmotif_merge:
    input:
        xmls=expand(join(WORK, "qmotif/{sample}/{sample}.xml"), sample=SAMPLES),
        wgs_metrics=WGS_METRICS_GRCH37
    output:
        tsv=join(WORK, "qmotif_merge.tsv")
    params:
        job_name="qmotif_merge", log="qmotif_merge.log"
    run:
        import xml.etree.ElementTree as ET
        from xml.dom import minidom

        coverage = get_coverage(input.wgs_metrics)

        with open(output.tsv, "w") as f_out:
            print("\t".join(["Sample", "Num_reads/coverage"]), file=f_out)

            for fname in input.xmls:
                subject_id = splitext(basename(fname))[0]

                line = [subject_id]

                with open(fname) as f_in:
                    xml = minidom.parse(f_in)

                    counts = xml.getElementsByTagName("rawIncludes")[0]

                    if counts.hasAttribute("count"):
                        length = float(counts.attributes["count"].value)
                    else:
                        length = 0

                    length_coverage = str(length / coverage[subject_id])

                    print("\t".join([subject_id, length_coverage]), file=f_out)

#=====================================================================================================
# motif_counter
#=====================================================================================================

rule motif_counter:
    input:
        bam=join(INPUT_DIR_BAM["grch38"], "{sample}.bam")
    output:
        join(WORK, "motif_counter/{sample}/out.txt")
    params:
        job_name="motif_counter-{sample}", log="motif_counter/{sample}.log", cpus="1", mem="5G", time="4-0",
        output_prefix=join(WORK, "motif_counter/{sample}")
    shell:
        """
        mkdir -p {params.output_prefix}
        cd {params.output_prefix}

        echo -e "TTAGGG\n9" | \
            {SINGULARITY} bash {BASE}/metadata/motif_counter.sh \
                -p -v \
                -i {input.bam} \
                -s -q 0 -Q 0 \
                -o out
        """

rule motif_counter_calc:
    input:
        tsvs=expand(join(WORK, "motif_counter/{sample}/out.txt"), sample=SAMPLES),
        wgs_metrics=WGS_METRICS_GRCH38
    output:
        tsv=join(WORK, "motif_counter.tsv")
    params:
        job_name="motif_counter_calc", log="motif_counter_calc.log"
    run:
        coverage = get_coverage(input.wgs_metrics)

        with open(output.tsv, "w") as f_out:
            print("\t".join(["Sample", "Num_reads/coverage"]), file=f_out)

            for fname in input.tsvs:
                with open(fname) as f_in:
                    next(f_in)

                    for line in f_in:
                        fields = line.split("\t")
                        subject_id, reads = fields[0], fields[4]
                        reads_coverage = str(float(reads) / coverage[subject_id])
                        print("\t".join([subject_id, reads_coverage]), file=f_out)

#=====================================================================================================
# combine methods
#=====================================================================================================

rule combine_methods:
    input:
        methods=[join(WORK, m) for m in COMBINE_METHODS]
    output:
        combined_methods=join(WORK, "combined_methods.tsv")
    run:
        from decimal import Decimal

        methods = [splitext(basename(m))[0] for m in input.methods]
        combined = {}

        for n, fname in enumerate(input.methods):
            method_name = methods[n]

            # <method_name>: [<0:separator>, <1:position_in_data>, <2:transformation_function>, <3:short_name>]
            method = {
                "computel_merge": ["\t", 0, lambda x: str(Decimal(x) / 1000), "computel"],
                "motif_counter": ["\t", 0, lambda x: x, "motif_counter"],
                "qmotif_merge": ["\t", 0, lambda x: x, "qmotif"],
                "telomerecat_telbam2length": [",", -1, lambda x: str(Decimal(x) / 1000), "telomerecat"],
                "telseq_merge": ["\t", 3, lambda x: x, "telseq"]
            }

            with open(fname) as f:
                next(f)
                for line in f:
                    sample, *data = line.rstrip("\n").split(method[method_name][0])
                    value = data[method[method_name][1]]
                    value = method[method_name][2](value)

                    combined.setdefault(sample, []).append(value)

        with open(output.combined_methods, "w") as f:
            print("\t".join(["sample"] + [method[m][3] for m in methods]), file=f)

            for sample in sorted(combined.keys()):
                print("\t".join([sample] + [v for v in combined[sample]]), file=f)

#=====================================================================================================
# correlation
#=====================================================================================================

rule correlation:
    input:
        combined_methods=join(WORK, "combined_methods.tsv")
    output:
        pearson_pdf=join(WORK, "correlation/pearson.pdf"),
        spearman_pdf=join(WORK, "correlation/spearman.pdf"),
        pearson_txt=join(WORK, "correlation/pearson.txt"),
        spearman_txt=join(WORK, "correlation/spearman.txt"),
        telomere_txt=join(WORK, "correlation/telomere.txt")
    params:
        job_name="correlation", log="correlation.log", output_prefix=join(WORK, "correlation/")
    shell:
        """
        {RSCRIPT} {BASE}/scripts/telomere_correlation.R {input.combined_methods} {params.output_prefix}
        """

# vim: ft=snakemake ts=4 sw=4 expandtab
