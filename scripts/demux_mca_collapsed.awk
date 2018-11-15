#!/usr/bin/awk -f
BEGIN {
    FS  = "\t";
    OFS = "\n";
    if (length(outdir) == 0 || length(prefix) == 0) {
	print "demux_mca.awk: Variables 'outdir' and 'prefix' must be defined!" > "/dev/stderr";
	exit 1;
    }

    outBase = outdir "/" prefix ".";
    print "[INFO] Initiating processing..." > "/dev/stderr";
}
{
    bx  = $1;
    umi = $2;
    
    if ( bx != prevbx ) {
	close(umiOut);
	close(fastqOut);
	umiOut   = outBase bx ".umi";
	fastqOut = outBase bx ".fastq";
	prevBX = bx;
    }

    print umi >> umiOut;
    print $3, $4, $5, $6 >> fastqOut;
    
    if (NR%100000 == 0) {
	printf("[INFO] %d reads processed\n", NR) > "/dev/stderr";
    }
}
END {
    printf("[INFO] %d reads processed\n", NR) > "/dev/stderr";
}
