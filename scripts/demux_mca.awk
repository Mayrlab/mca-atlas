#!/usr/bin/awk -f
BEGIN {
    if (length(outdir) == 0 || length(prefix) == 0) {
	print "demux_mca.awk: Variables 'outdir' and 'prefix' must be defined!" > "/dev/stderr";
	exit 1;
    }
    print "[INFO] Initiating processing..." > "/dev/stderr";
}
{
    if (NR%4 == 1) {
	match($1, /.*_([ACGT]{18})_([ACGTN]{6}).*/, bx);
	print bx[2] >> outdir"/"prefix"."bx[1]".umi";
    }
    print >> outdir"/"prefix"."bx[1]".fastq";
    
    if (NR%40000 == 0) {
	printf("[INFO] %d reads processed\n", NR/4) > "/dev/stderr";
    }
}
END {
    printf("[INFO] %d reads processed\n", NR/4) > "/dev/stderr";
}
