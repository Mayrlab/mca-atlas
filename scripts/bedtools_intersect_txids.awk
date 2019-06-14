#!/usr/bin/gawk -f
#' Extract transcript_name -> transcript_name map of overlaps from bedtools intersect
BEGIN {
    OFS = "\t"
}
{
    match($0, /transcript_id "([^"]+\.[0-9]+)".*transcript_id "([^"]+\.[0-9]+)"/, txids);
    if (txids[1] != txids[2]) {
        print txids[1], txids[2];
    } 
}
