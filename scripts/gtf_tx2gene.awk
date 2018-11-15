#!/usr/bin/gawk -f
#' Extract TX -> GENE map from GTF, including chromosome
BEGIN {
    OFS = "\t"
}
{
    if ($3 ~ /transcript/) {
	match($0, /gene_id "(ENSMUSG[0-9]+\.[0-9]+)";.*transcript_id "([^"]+\.[0-9]+)"/, ids);
	match(ids[2], /(.*)\.[0-9]+$/, symbol);
	print ids[2], ids[1], symbol[1], $1;
    }
}
