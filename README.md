## 3' UTR cleavage site identification from the Mouse Cell Atlas
This reposistory hosts a Snakemake pipeline used to call 3' UTR cleavage sites *de novo* from [the Mouse Cell Atlas dataset](http://bis.zju.edu.cn/MCA/). The final output of the pipeline is the Mouse UTRome used in [the scUTRquant pipeline](https://github.com/Mayrlab/scUTRquant) for quantifying 3' UTR isoforms from scRNA-seq data.

Please see [our accompanying *bioRxiv* manuscript](https://doi.org/10.1101/2021.11.22.469635).

## Requirements
The original pipeline was run with:

- Snakemake 5.10
- kallisto 0.44.0

Note that on LSF system, we additionally used the `scripts/snakemake_bsub.py` to manage cluster scaling. Rule resources may need to be adjusted for other systems.

## Running
After cloning this repository, a kallisto index for a UTRome can be generated with:

```bash
snakemake data/kallisto/adult.utrome.e3.t200.f0.9999.w500.kdx
```

Additional parameters are encoded in the filename:

 - **epsilon (e):** distance for initial merging of cleavage sites
 - **threshold (t):** minimum number of supporting reads
 - **non-internal priming likelihood (f):** minimum posterior likelihood for not being an internal priming site
 - **truncation width (w):** maximum transcript distance to cleavage site
