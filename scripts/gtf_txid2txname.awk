#!/usr/bin/gawk -f
#' Extract transcript_id -> transcript_name map from GTF, including chromosome
BEGIN {
    OFS = "\t"
}
{
    if ($3 ~ /transcript/) {
	match($0, /transcript_id "([^"]+\.[0-9]+)";.*transcript_name "([^"]+\.[0-9]+(-UTR(-|\+)[0-9]+|))"/, ids);
	match(ids[1], /(.*)\.[0-9]+$/, symbol);
	if (ids[2] ~ /UTR-/) {
	    utr = "upstream"
	} else if (ids[2] ~ /UTR\+/) {
	    utr = "downstream"
	} else {
	    utr = "GENCODE"
	}
	print ids[1], ids[2], symbol[1], $1, utr;
    }
}
