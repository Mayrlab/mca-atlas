configfile: "config.yaml"

wildcard_constraints:
    srr="SRR\d+"

shell.prefix("source ~/.bashrc; ")
    
import pandas as pd
import os

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)

print("Loading annotation data...")
annotations = pd.read_csv(config["annotationFile"])

metadata = pd.read_table(config["metadataFile"], usecols=["Run", "Sample_Name", "tissue"])
metadata["tissue_srr"] = metadata.tissue.map(str) + "." + metadata.Run

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
adultSRRs = adultAnnotations.Run.unique()

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

workingSRRs = adultSRRs # config["immuneSRRs"]

rule all:
    input:
        expand("qc/raw/{srr}_{read}_fastqc.html", srr = workingSRRs, read = [1,2]),
        expand("data/fastq/assembled/{srr}.assembled.fastq.gz", srr = workingSRRs),
        expand("qc/trimmed/{srr}.{readtype}.clean_fastqc.html", srr = workingSRRs, readtype = ["assembled", "unassembled_1", "unassembled_2"]),
        expand("output/read_distribution/{tissue_srr}.{readtype}.png",
               tissue_srr = metadata.loc[metadata["Run"].isin(workingSRRs), "tissue_srr"], readtype = ["assembled", "unassembled"]),
        expand("output/read_distribution/MouseCellAtlas.AdultTissues.{readtype}.png", readtype = ["assembled", "unassembled"]),
        expand("data/coverage/adult.assembled.{strand}.txt.gz", strand=['positive', 'negative'])
               

rule download:
    output:
        temp("data/sra/{srr}.sra")
    params:
        ftpbase = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR",
        filepath = lambda wcs: "{srr6}/{srr}/{srr}.sra".format(srr6 = wcs.srr[0:6], srr = wcs.srr),
        bitrate = "1000m",
        sshkey = config["ascpKeyFile"]
    shell:
        """
        ascp -i {params.sshkey} -k 2 -T -l{params.bitrate} {params.ftpbase}/{params.filepath} {output}
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

rule fastqc_clean:
    input:
        "data/fastq/trimmed/{srr}.{readtype}.clean.fastq.gz"
    output:
        "qc/trimmed/{srr}.{readtype}.clean_fastqc.html",
        "qc/trimmed/{srr}.{readtype}.clean_fastqc.zip"
    params:
        adapters = config["adaptersFile"],
        tmp_dir = config["tmp_dir"]
    shell:
        """
        fastqc -a {params.adapters} -d {params.tmp_dir} -o qc/trimmed {input}
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
    threads: 24
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
        mem = 16,
        walltime = 24
    shell:
        """
        umi_tools whitelist --set-cell-number={params.cell_num} --method=umis \\
        --extract-method=regex --plot-prefix={params.prefix} \\
        --bc-pattern='(?P<cell_1>.{{6}})(?P<discard_1>CGACTCACTACAGGG){{s<=1}}(?P<cell_2>.{{6}})(?P<discard_2>TCGGTGACACGATCG){{s<=1}}(?P<cell_3>.{{6}})(?P<umi_1>.{{6}})(T{{12}}){{s<=2}}.*' \\
        --stdin={input} --stdout={output.raw} --log={log}
        awk '$3 >= 500 {{ print $0 }}' {output.raw} > {output.filtered}
        """

rule umitools_extract_assembled:
    input:
        bx = "data/barcodes/{srr}.raw.whitelist.umi500.txt",
        fq = "data/fastq/assembled/{srr}.assembled.fastq.gz"
    output:
        "data/fastq/extracted/{srr}.assembled.bx.fastq.gz"
    log:
        "logs/umi_tools/{srr}.assembled.extract.log"
    resources:
        mem = 16,
        walltime = 24
    shell:
        """
        umi_tools extract --filter-cell-barcode --error-correct-cell \\
        --extract-method=regex \\
        --bc-pattern='(?P<cell_1>.{{6}})(?P<discard_1>CGACTCACTACAGGG){{s<=1}}(?P<cell_2>.{{6}})(?P<discard_2>TCGGTGACACGATCG){{s<=1}}(?P<cell_3>.{{6}})(?P<umi_1>.{{6}})(T{{12}}){{s<=2}}.*' \\
        --whitelist={input.bx} --log={log} \\
        --stdin={input.fq} --stdout={output}
        """

rule umitools_extract_unassembled:
    input:
        bx = "data/barcodes/{srr}.raw.whitelist.umi500.txt",
        r1 = "data/fastq/unassembled/{srr}.unassembled_1.fastq.gz",
        r2 = "data/fastq/unassembled/{srr}.unassembled_2.fastq.gz"
    output:
        r1 = "data/fastq/extracted/{srr}.unassembled_1.bx.fastq.gz",
        r2 = "data/fastq/extracted/{srr}.unassembled_2.bx.fastq.gz"
    log:
        "logs/umi_tools/{srr}.unassembled.extract.log"
    resources:
        mem = 16,
        walltime = 24
    shell:
        """
        umi_tools extract --filter-cell-barcode --error-correct-cell \\
        --extract-method=regex \\
        --bc-pattern='(?P<cell_1>.{{6}})(?P<discard_1>CGACTCACTACAGGG){{s<=1}}(?P<cell_2>.{{6}})(?P<discard_2>TCGGTGACACGATCG){{s<=1}}(?P<cell_3>.{{6}})(?P<umi_1>.{{6}})(T{{12}}){{s<=2}}.*' \\
        --whitelist={input.bx} --log={log} \\
        --stdin={input.r1} --read2-in={input.r2} \\
        --stdout={output.r1} --read2-out={output.r2}
        """

rule cutadapt_polyT_SE:
    input:
        "data/fastq/extracted/{srr}.assembled.bx.fastq.gz"
    output:
        "data/fastq/trimmed/{srr}.assembled.clean.fastq.gz"
    params:
        min_length = config["minReadLength"]
    threads: 8
    shell:
        """
        cutadapt --cores={threads} --front='T{{100}}' \\
        --minimum-length={params.min_length} --length-tag='length=' \\
        --output={output} {input}
        """

rule cutadapt_polyT_PE:
    input:
        r1 = "data/fastq/extracted/{srr}.unassembled_1.bx.fastq.gz",
        r2 = "data/fastq/extracted/{srr}.unassembled_2.bx.fastq.gz"
    output:
        r1 = "data/fastq/trimmed/{srr}.unassembled_1.clean.fastq.gz",
        r2 = "data/fastq/trimmed/{srr}.unassembled_2.clean.fastq.gz",
    params:
        min_length = config["minReadLength"]
    threads: 8
    shell:
        """
        cutadapt --cores={threads} --front='T{{100}}' -A='A{{100}}' \\
        --minimum-length={params.min_length} --pair-filter=both --length-tag='length=' \\
        --output={output.r1} --paired-output={output.r2} {input.r1} {input.r2}
        """

rule read_distribution:
    input:
        bam = "data/bam/{srr}.{readtype}.bam",
        bed = config["gencodeBED"]
    output:
        "qc/aligned/{srr}.{readtype}.read_distribution.txt"
    shell:
        """
        source activate py27
        read_distribution.py -i {input.bam} -r {input.bed} &> {output}
        """

rule read_dist_plot:
    input:
        "qc/aligned/{srr}.{readtype}.read_distribution.txt"
    output:
        "output/read_distribution/{tissue}.{srr}.{readtype}.png"
    params:
        label = lambda wcs: '"Mouse Cell Atlas {tissue} Sample ({srr}, {readtype})"'.format(tissue = wcs.tissue, srr = wcs.srr, readtype = wcs.readtype)
    shell:
        """
        ./scripts/read_dist_plot.R {params.label} {input} {output}
        """

rule hisat2_PE:
    input:
        r1 = "data/fastq/trimmed/{srr}.unassembled_1.clean.fastq.gz",
        r2 = "data/fastq/trimmed/{srr}.unassembled_2.clean.fastq.gz"
    output:
        bam = "data/bam/{srr}.unassembled.bam",
        bai = "data/bam/{srr}.unassembled.bam.bai"
    params:
        sam = config["tmp_dir"] + "/{srr}.unassembled.sam",
        tmp_dir = config["tmp_dir"],
        idx = config['hisatIndex']
    threads: 16
    log: "logs/hisat2/{srr}.unassembled.log"
    shell:
        """
        hisat2 -p {threads} -x {params.idx} -1 {input.r1} -2 {input.r2} -S {params.sam} --new-summary --summary-file {log}
        samtools sort -@ {threads} -m 6G -T {params.tmp_dir}/ -o {output.bam} {params.sam}
        samtools index -@ {threads} {output.bam}
        rm -f {params.sam}
        """

rule hisat2_SE:
    input:
        "data/fastq/trimmed/{srr}.assembled.clean.fastq.gz"
    output:
        bam = "data/bam/{srr}.assembled.bam",
        bai = "data/bam/{srr}.assembled.bam.bai"
    params:
        tmp_dir = config["tmp_dir"],
        idx = config['hisatIndex'],
        sam = config["tmp_dir"] + "/{srr}.assembled.sam"
    threads: 16
    log: "logs/hisat2/{srr}.assembled.log"
    shell:
        """
        hisat2 -p {threads} -x {params.idx} -U {input} -S {params.sam} --new-summary --summary-file {log}
        samtools sort -@ {threads} -m 6G -T {params.tmp_dir}/ -o {output.bam} {params.sam}
        samtools index -@ {threads} {output.bam}
        rm -f {params.sam}
        """

rule read_dist_plot_adult:
   input:
       expand("qc/aligned/{srr}.{readtype}.read_distribution.txt", srr = adultSRRs, readtype = ['assembled', 'unassembled'])
   output:
       assembled="output/read_distribution/MouseCellAtlas.AdultTissues.assembled.png",
       unassembled="output/read_distribution/MouseCellAtlas.AdultTissues.unassembled.png"
   params:
       metadataFile = config["metadataFile"],
       rx_assembled = "'SRR\\d+.assembled.read_distribution.txt$'",
       rx_unassembled = "'SRR\\d+.unassembled.read_distribution.txt$'"
   shell:
       """
       ./scripts/read_dist_plot_all.R {params.metadataFile} qc/aligned {params.rx_assembled} {output.assembled}
       ./scripts/read_dist_plot_all.R {params.metadataFile} qc/aligned {params.rx_unassembled} {output.unassembled}
       """

rule cleavage_coverage:
    input:
        "data/bam/{srr}.assembled.bam"
    output:
        neg = "data/coverage/{srr}.assembled.negative.txt.gz",
        pos = "data/coverage/{srr}.assembled.positive.txt.gz"
    shell:
        """
        bedtools genomecov -dz -5 -strand '-' -ibam {input} | gzip > {output.neg}
        bedtools genomecov -dz -5 -strand '+' -ibam {input} | gzip > {output.pos}
        """

rule sum_cov_adult:
    input:
        neg = expand("data/coverage/{srr}.assembled.negative.txt.gz", srr=adultSRRs),
        pos = expand("data/coverage/{srr}.assembled.positive.txt.gz", srr=adultSRRs)
    output:
        neg = "data/coverage/adult.assembled.negative.txt.gz",
        pos = "data/coverage/adult.assembled.positive.txt.gz"
    threads: 4
    shell:
        """
        zcat {input.neg} | datamash -s -g 1,2 sum 3 | sort -k 1,1 -k 2,2n | gzip > {output.neg}
        zcat {input.pos} | datamash -s -g 1,2 sum 3 | sort -k 1,1 -k 2,2n | gzip > {output.pos}
        """

rule merge_cov_adult:
    input:
        "data/coverage/adult.assembled.{strand}.txt.gz"
    output:
        "data/coverage/adult.assembled.{strand}.{epsilon,\d+}.txt.gz"
    threads: 16
    resources:
        walltime=24
    shell:
        """
        scripts/merge_genomecov.py -p {threads} -e {wildcards.epsilon} {input} {output}
        """

rule cov_to_bed:
    input:
        neg = "data/coverage/adult.assembled.negative.{epsilon}.txt.gz",
        pos = "data/coverage/adult.assembled.positive.{epsilon}.txt.gz"
    output:
        gz = "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed.gz",
        tbi = "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed.gz.tbi"
    params:
        raw = "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+"
    shell:
        """
        zcat {input.neg} | awk -v OFS='\\t' '$3>={wildcards.threshold}{{print $1, $2, $2, $1 \":\" $2 \":+\", $3, \"+\"}}' > {params.raw}
        zcat {input.pos} | awk -v OFS='\\t' '$3>={wildcards.threshold}{{print $1, $2, $2, $1 \":\" $2 \":-\", $3, \"-\"}}' >> {params.raw}
        sort -k1,1 -k2,2n {params.raw} | bgzip > {output.gz}
        tabix -0 -p bed {output.gz}
        rm -f {params.raw}
        """

rule cleanUpdTSeq_adult:
    input:
        "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed.gz"
    output:
        res = "qc/cleavage-sites/adult.classified.e{epsilon}.t{threshold}.tsv.gz",
        hist = "qc/cleavage-sites/adult.posteriors.e{epsilon}.t{threshold}.png"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+"
    resources:
        walltime=24,
        mem = 24
    shell:
        """
        scripts/classify_cleavage_sites.R {input} {output.res} {output.hist}
        """

rule intersect_polyASite:
    input:
        "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed.gz"
    output:
        "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.pA{score}.bed.gz"
    params:
        pasFile=config["polyASiteBED"],
        window=12
    shell:
        """
        scripts/compare_polyASite.R {params.pasFile} {wildcards.score} {params.window} {input} {output}
        """

rule export_cleavage_sites:
    input:
        full = "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.bed.gz",
        clean = "data/bed/cleavage-sites/adult.cleavage.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        pas = config["polyASiteBED"],
        annot = config["gencodeGFF"]
    output:
        "data/bed/cleavage-sites/adult.validated.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "data/bed/cleavage-sites/adult.supported.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "data/bed/cleavage-sites/adult.likely.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "data/bed/cleavage-sites/adult.unlikely.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "data/bed/cleavage-sites/adult.utr3.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "data/bed/cleavage-sites/adult.extutr3.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        "qc/cleavage-sites/adult.annotations.e{epsilon}.t{threshold}.f{likelihood}.png",
        "qc/cleavage-sites/adult.site_scores.e{epsilon}.t{threshold}.f{likelihood}.png"
    threads: 1
    params:
        prefix = "adult",
        suffix = lambda wcs: "e%s.t%s.f%s" % (wcs.epsilon, wcs.threshold, wcs.likelihood),
        extUTR3 = 5000,
        extUTR5 = 1000
    shell:
        """
        scripts/export_cleavage_sites.R {threads} {input.full} {input.clean} {input.annot} {input.pas} {params.extUTR5} {params.extUTR3} {params.prefix} {params.suffix}
        """

rule augment_transcriptome:
    input:
        utr3 = "data/bed/cleavage-sites/adult.utr3.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        ext3 = "data/bed/cleavage-sites/adult.extutr3.e{epsilon}.t{threshold}.f{likelihood}.bed.gz",
        annot = config["gencodeGFF"]
    output:
        "data/gff/adult.txs.utr3.e{epsilon}.t{threshold}.f{likelihood}.gff3",
        "data/gff/adult.txs.utr3.e{epsilon}.t{threshold}.f{likelihood}.gtf",
        "data/gff/adult.txs.extutr3.e{epsilon}.t{threshold}.f{likelihood}.gff3",
        "data/gff/adult.txs.extutr3.e{epsilon}.t{threshold}.f{likelihood}.gtf"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+",
        likelihood = "0.\d+"
    threads: 24
    resources:
        mem = 8,
        walltime = 24
    params:
        prefix = "adult",
        suffix = lambda wcs: "e%s.t%s.f%s" % (wcs.epsilon, wcs.threshold, wcs.likelihood),
        extUTR3 = 5000
    shell:
        """
        scripts/augment_transcriptome.R {threads} {input.utr3} {input.ext3} {input.annot} {params.extUTR3} {params.prefix} {params.suffix}
        """

rule export_utrome:
    input:
        utr3 = "data/gff/adult.txs.utr3.e{epsilon}.t{threshold}.f{likelihood}.gff3",
        ext3 = "data/gff/adult.txs.extutr3.e{epsilon}.t{threshold}.f{likelihood}.gff3",
        annot = config["gencodeSortedGFF"]
    output:
        gtf = "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.gtf",
        fa = "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.fasta"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+",
        likelihood = "0.\d+",
        width = "\d+"
    threads: 4
    params:
        prefix = "adult",
        suffix = lambda wcs: "e%s.t%s.f%s.w%s" % (wcs.epsilon, wcs.threshold, wcs.likelihood, wcs.width)
    shell:
        """
        scripts/export_utrome.R {threads} {input.utr3} {input.ext3} {input.annot} {wildcards.width} {params.prefix} {params.suffix}
        """

rule sort_utrome:
    input:
        "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.gtf"
    output:
        gtf = "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.gtf.gz",
        idx = "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.gtf.gz.tbi"
    threads: 4
    resources:
        mem = 6
    shell:
        """
        cat {input} | awk '!( $0 ~ /^#/ )' | sort --parallel={threads} -S4G -k1,1 -k4,4n | bgzip -c > {output.gtf}
        tabix {output.gtf}
        """

rule plot_utrome_stats:
    input:
        "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.gtf.gz"
    output:
        "qc/utrome/adult.utrome.txsPerGene.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.exonsPerTx.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.lengths.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.dists.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.distsByChr.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.sitesPerGene.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.sitesPerGene.log.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.txsPerSite.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png",
        "qc/utrome/adult.utrome.txsPerSite.log.e{epsilon}.t{threshold}.f{likelihood}.w{width}.png"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+",
        likelihood = "0.\d+",
        width = "\d+"
    threads: 1
    params:
        prefix = "adult",
        suffix = lambda wcs: "e%s.t%s.f%s.w%s" % (wcs.epsilon, wcs.threshold, wcs.likelihood, wcs.width)
    shell:
        """
        scripts/plot_utrome_stats.R {threads} {input} {params.prefix} {params.suffix}
        """

        
rule kallisto_index:
    input:
        "data/gff/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.fasta"
    output:
        "data/kallisto/adult.utrome.e{epsilon}.t{threshold}.f{likelihood}.w{width}.kdx"
    resources:
        mem = 16
    shell:
        """
        kallisto index -i {output} {input}
        """

rule bam_mv_bxs:
    input:
        "data/bam/{srr}.{readtype}.bam"
    output:
        bam = "data/bam/{srr}.{readtype}.bxs.bam",
        idx = "data/bam/{srr}.{readtype}.bxs.bam.bai"
    shell:
        """
        samtools view -h {input} | awk '{ if($0 ~ "^@") {print $0} else { split($1, read_name, "_"); $1 = read_name[1]; print $0"\tCB:Z:"read_name[2]"\tRX:Z:"read_name[3]}}' | samtools view -b > {output.bam}
        samtools index {output.bam}
        """

##rule demux_map:
##    input:
##        "data/fastq/trimmed/{srr}.{readtype}.clean.fastq.gz"
##    output:
##        bam = "data/
