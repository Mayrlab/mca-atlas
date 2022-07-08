## 3' UTR cleavage site identification from the Mouse Cell Atlas
This reposistory hosts the Snakemake pipeline used to call 3' UTR cleavage sites *de novo* from the Mouse Cell Atlas dataset. The final output of the pipeline is the Mouse UTRome used in the scUTRquant pipeline for quantifying 3' UTR isoforms from scRNA-seq data.

## Requirements
The original pipeline was run with:

- Snakemake 5.10
- kallisto 0.44.0

## Running
After cloning, a kallisto index for a UTRome can be generated with:

```bash
snakemake data/kallisto/adult.utrome.e3.t200.f0.9999.w500.kdx
```

Additional parameters are encoded in the filename:

 - **epsilon (e):** distance for initial merging of cleavage sites
 - **threshold (t):** minimum number of supporting reads
 - **non-internal priming likelihood (f):** minimum posterior likelihood of being an internal priming site
 - **truncation width (w):** maximum transcript distance to cleavage site
