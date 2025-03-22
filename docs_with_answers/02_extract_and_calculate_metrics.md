# Extraction and calculation answers

The answers are based on the files provided within the repo. If you have files customised yourself, the answers may very depending on the question. In some cases, there are more ways to obtain the answer, the answers provided are one way to do this. Sometimes the answers outlined below are not directly an answer to the questions in the exercise, but you will be able to manually count or infer the answers from them.


**Task: Generate individual metrics**

Extract GQ and DP information can be done either by generating one file per metric, or two-metrics-in-one file. We recommend the latter to remain more storage efficient.

```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
OUTPUT_FILE="aux/raw_sample_variant_gq_dp.tsv"

bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GQ\t%DP\n]\n' ${INPUT_FILE} > ${OUTPUT_FILE}
```
We do not provide the `aux/raw_sample_variant_gq_dp.tsv` file on our repo to save on storage (~2 Gb).

**Task: Calculate medianGQ and medianDP**

Calculating the medianGQ and medianDP value per variant was done through the script in `bin/calculate_median_gq_dp.R`. The output of this script is: `aux/gq_dp_annotation_file.txt` (~1.6 Mb).

1. What are the highest and lowest raw GQ and DP values you observe?
	- Min raw GQ: `.` or `0`
	- Max raw GQ: `99`
	- Min raw DP: `.` or `0`
	- Max raw DP: `99`
2. What are the highest and lowest median GQ and DP values you observe?
	- Min median GQ: `0`
	- Max median GQ: `99`
	- Min median DP: `0`
	- Max median DP: `56`

