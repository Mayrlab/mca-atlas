configfile: "config.yaml"

wildcard_constraints:
    srr="SRR\d+"

import pandas as pd
import os

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)

print("Loading annotation data...")
annotations = pd.read_csv(config["annotationFile"])

metadata = pd.read_table(config["metadataFile"], usecols=["Run", "Sample_Name", "tissue"])

sample_batch = pd.read_table(config["sampleBatchFile"])
metadata = pd.merge(metadata, sample_batch, on = "Sample_Name")

annotations = pd.merge(annotations, metadata, on = "Batch")

print("Loaded annotations for %d cells from %d libraries and %d tissues, comprising %d clusters." % (
    len(annotations), len(annotations.Batch.unique()), len(annotations.Tissue.unique()), len(annotations.ClusterID.unique())))

adults = [
    '8-10 week-old male mice',
    ' 6-10 week-old female mice',
    '8-11 week-old male mice']

adultAnnotations = annotations.loc[annotations["Mouse-Sex-Age"].isin(adults)]

print("This included:\n\t%d adult cells from %d libraries and %d tissues, comprising %d clusters;" % (
    len(adultAnnotations), len(adultAnnotations.Batch.unique()), len(adultAnnotations.Tissue.unique()), len(adultAnnotations.ClusterID.unique())))

neonatalAnnotations = annotations.loc[annotations["Mouse-Sex-Age"] == " 1 day-old pups "]

print("\t%d neonatal cells from %d libraries and %d tissues, comprising %d clusters." % (
    len(neonatalAnnotations), len(neonatalAnnotations.Batch.unique()), len(neonatalAnnotations.Tissue.unique()), len(neonatalAnnotations.ClusterID.unique())))

embryonicAnnotations = annotations.loc[annotations["Mouse-Sex-Age"] == "E14.5 embryos"]

print("\t%d embryonic cells from %d libraries and %d tissues, comprising %d clusters." % (
    len(embryonicAnnotations), len(embryonicAnnotations.Batch.unique()), len(embryonicAnnotations.Tissue.unique()), len(embryonicAnnotations.ClusterID.unique())))

stemCells = [
    'E14 cell line ',
    ' Mesenchymal stem cell line',
    'Trophoblast stem cell line']

stemCellAnnotations = annotations.loc[annotations["Mouse-Sex-Age"].isin(stemCells)]

print("\t%d stem cells from %d libraries and %d tissues, comprising %d clusters." % (
    len(stemCellAnnotations), len(stemCellAnnotations.Batch.unique()), len(stemCellAnnotations.Tissue.unique()), len(stemCellAnnotations.ClusterID.unique())))

rule all:
    input:
        expand("qc/raw/{srr}_{read}_fastqc.html", srr = config["immuneSRRs"], read = [1,2]),
        expand("data/fastq/assembled/{srr}.assembled.fastq.gz", srr = config["immuneSRRs"])

rule download:
    output:
        temp("data/sra/{srr}.sra")
    params:
        ftpbase = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR",
        filepath = lambda wcs: "{srr6}/{srr}/{srr}.sra".format(srr6 = wcs.srr[0:5], srr = wcs.srr),
        bitrate = "1000m",
        sshkey = config["ascpKeyFile"]
    shell:
        """
        ascp -i {params.sshkey} -k 2 -T -l{params.bitrate} {params.ftpbase}/{params.filepath}
        """

rule fastqdump_sra:
    input:
        "data/sra/{srr}.sra"
    output:
        "data/fastq/raw/{srr}_1.fastq.gz",
        "data/fastq/raw/{srr}_2.fastq.gz"
    params:
        tmp_dir = config["tmp_dir"]
    threads: 8
    shell:
        """
        parallel-fastq-dump -s {input} -t {threads} -O data/fastq/raw --tmpdir {params.tmp_dir} --split-files --gzip
        """

rule fastqc_raw:
    input:
        "data/fastq/raw/{srr}_{read}.fastq.gz"
    output:
        "qc/raw/{srr}_{read}_fastqc.html",
        "qc/raw/{srr}_{read}_fastqc.zip"
    params:
        adapters = config["adaptersFile"],
        tmp_dir = config["tmp_dir"]
    shell:
        """
        fastqc -a {params.adapters} -d {params.tmp_dir} -o qc/raw {input}
        """

rule pear_merge:
    input:
        r1 = "data/fastq/raw/{srr}_1.fastq.gz",
        r2 = "data/fastq/raw/{srr}_2.fastq.gz"
    output:
        r12 = "data/fastq/assembled/{srr}.assembled.fastq.gz",
        r1 = "data/fastq/unassembled/{srr}.unassembled_1.fastq.gz",
        r2 = "data/fastq/unassembled/{srr}.unassembled_2.fastq.gz",
        r0 = "data/fastq/discarded/{srr}.discarded.fastq.gz"
    params:
        tmp_dir = config["tmp_dir"] + "/pear",
        total_mem = 32,
        min_length = 54 + 21
    threads: 8
    log:
        "logs/pear/{srr}.log"
    shell:
        """
        mkdir -p {params.tmp_dir}
        pear -j {threads} -n {params.min_length} -p 0.0001 -f {input.r1} -r {input.r2} -o {params.tmp_dir}/{wildcards.srr} 2> {log}
        pigz -p {threads} -c {params.tmp_dir}/{wildcards.srr}.assembled.fastq > {output.r12}
        pigz -p {threads} -c {params.tmp_dir}/{wildcards.srr}.unassembled.forward.fastq > {output.r1}
        pigz -p {threads} -c {params.tmp_dir}/{wildcards.srr}.unassembled.reverse.fastq > {output.r2}
        pigz -p {threads} -c {params.tmp_dir}/{wildcards.srr}.discarded.fastq > {output.r0}
        rm {params.tmp_dir}/{wildcards.srr}.*.fastq
        """

# rule umitools_whitelist:
#     input:
#         r12 = "data/fastq/assembled/{srr}.assembled.fastq.gz",
#         r1 = "data/fastq/unassembled/{srr}.unassembled_1.fastq.gz",
#         r2 = "data/fastq/unassembled/{srr}.unassembled_2.fastq.gz"
#     output:
#         r12 = "data/barcodes/{srr}.assembled.whitelist.txt",
#         r1 = "data/barcodes/{srr}.unassembled.whitelist.txt"
#     threads: 4
#     params:
#         bc_regex = config["barcodeRegex"]
#     log:
#         r12 = "logs/umi_tools/{ssr}.assembled.whitelist.log",
#         r1 = "logs/umi_tools/{ssr}.unassembled.whitelist.log"
#     shell:
#         """
#         umi_tools whitelist -I {input.r12} -S {output.r12} -L {log.r12} --extract-method='regex' --bc-pattern='{params.bc_regex}'
#         umi_tools whitelist -I {input.r1} -S {output.r1} -L {log.r1} --extract-method='regex' --bc-pattern='{params.bc_regex}'
#         """

rule umitools_whitelist_raw:
    input:
        "data/fastq/raw/{srr}_1.fastq.gz"
    output:
        raw = "data/barcodes/{srr}.raw.whitelist.c20K.txt",
        filtered = "data/barcodes/{srr}.raw.whitelist.umi500.txt"
    params:
        #bc_regex = config["barcodeRegex"],
        cell_num = 20000,
        prefix = lambda wcs: "qc/barcodes/%s.raw.20K" % wcs.srr
    log:
        "logs/umi_tools/{srr}.raw.whitelist.c20K.log"
    resources:
        walltime = 24
    shell:
        """
        umi_tools whitelist --set-cell-number={params.cell_num} --method=umis \\
        --extract-method=regex --plot-prefix={params.prefix} \\
        --bc-pattern='(?P<cell_1>.{{6}})(?P<discard_1>CGACTCACTACAGGG){{s<=1}}(?P<cell_2>.{{6}})(?P<discard_2>TCGGTGACACGATCG){{s<=1}}(?P<cell_3>.{{6}})(?P<umi_1>.{{6}})(T{{12}}){{s<=2}}.*' \\
        --stdin={input} --stdout={output.raw} --log={log}
        awk '$3 >= 500 {{ print $0 }}' {output.raw} > {output.filtered}
        """
        
#rule cutadapt:
#    input:
#        "data/fastq/raw/{srr}_{read}.fastq.gz"
#    output:
#        "qc/trimmed/{srr}_{read}.trimmed.gz"

rule read_distribution:
    input:
        bam = "data/bam/{tissue}.{srr}.bam",
        bed = config["gencodeBED"]
    output:
        "qc/aligned/{tissue}.{srr}.read_distribution.txt"
    shell:
        """
        source activate py27
        read_distribution.py -i {input.bam} -r {input.bed} &> {output}
        """

rule read_dist_plot:
    input:
        "qc/aligned/{tissue}.{srr}.read_distribution.txt"
    output:
        "output/read_distribution/{tissue}.{srr}.png"
    params:
        label = lambda wcs: '"Mouse Cell Atlas {tissue} Sample ({srr})"'.format(tissue = wcs.tissue, srr = wcs.srr)
    shell:
        """
        ./scripts/read_dist_plot.R {params.label} {input} {output}
        """

#rule read_dist_plot_all:
#    input:
#        expand("qc/aligned/mca.{tissue}.{i}.read_distribution.txt", srr = list(metadata.Run))
#    output:
#        "output/read_distribution/TabulaMuris.10Xsamples.png"
#    params:
#        metadataFile = config["metadataFile"]
#    shell:
#        """
#        ./scripts/read_dist_plot_all.R {params.metadataFile} qc/aligned {output}
#        """
