#!/usr/bin/awk -f
BEGIN {
    FS="[_[:space:]]";
    OFS="\t";
    ORS="";
    print "[INFO] Collapsing FASTQ to single line format...\n" > "/dev/stderr";
}
{
    if (NR%4 == 1) {
	print $2 OFS $3 OFS $0;
    } else {
	print OFS $0;
    }
    if (NR%4 == 0) {
	print "\n";
    }
    
    if (NR%400000 == 0) {
	printf("[INFO] %d reads collapsed\n", NR/4) > "/dev/stderr";
    }
}
END {
    printf("[INFO] Collapsed %d reads\n", NR/4) > "/dev/stderr";
}
