# VCF data extraction and calculating QC metrics 

This section is directed at extracting information from a VCF, and calculating a new variant-level metric which will be integrated into the original VCF at a later stage.

## Exercise requirements

Along with knowledge of some basic linux commands, the following is needed for this exercise.
**Note:** You will also require data wrangling skills in `R` or `Python`, or any other method of choice which can calculate the median across a vector of values.

Required tools:
```
bcftools --version
```

Required (ready made) files:
```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
```

## Exercise

To ensure that we have confidence in our variant data, we want to know per variant what their quality roughly is across all samples. Sometimes variants that are located in repetitive, or high complexity regions, the quality levels may be low. Therefore, it helps researchers estimate whether to have confidence in a given variant when it is potentially marked as clinically significant. 

In this exercise, we will start preparing our quality assessment by first extracting the information we need, and then calculating variant-level metrics with the tool of your choice. Right now, genotype quality and read depth are only available as a value per sample. We will want to calculate the median GQ and median DP to reduce this information into a single value per variant, which will ultimately allow us to infer whether a variant is reliable.

Metrics of interest:
- `GQ`: Genotype quality, as a phred-scaled confidence on whether the genotype call is wrong is correct. While higher is better, it is often capped at 99.
- `DP`: Read depth at this position for the sample.

Tasks:
1. Extract for each variant and sample the genotype quality (`GQ`) and read depth (`DP`). 
2. Calculate per variant the median GQ and median DP.
	- What are the highest and lowest raw GQ and DP values you observe?
	- What are the highest and lowest median GQ and DP values you observe? 
3. Generate a `.tsv` file with: `CHROM,POS,REF,ALT,medianGQ,medianDP`


## Outcomes

You should now (roughly) know:
- Extract specific metrics data
- Build and prepare a file that can be used for annotation



