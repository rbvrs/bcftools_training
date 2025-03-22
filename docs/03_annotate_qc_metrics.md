# Annotating our VCF with QC metrics

This section is directed at annotating our VCF with QC-metrics that were calculated in the previous exercise.

## Exercise requirements

Along with knowledge of some basic linux commands, the following is needed for this exercise.

Required tools:
```
bcftools --version
bgzip --version
tabix --version
```

Required (ready made) files:
```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
ANNOTATION_FILE="aux/gq_dp_annotation_file.txt"
```
If you made your own annotation file, adjust the file paths accordingly.

## Exercise

To annotate a VCF with your own custom dataset, you will need three things: 1) Lines to add to the header, 2) a gzipped input annotation file, and 3) the target VCF. This exercise is aimed at learning how to annotate VCFs with your own custom data. 

1. Generate header lines for `medianGQ` and `medianDP`. Also think about whether these should these be considered `INFO` or `FORMAT` fields? And why?
2. Gzip and index the annotation file.
3. Annotate the VCF with `bcftools annotate`.
4. Confirm that your output VCF indeed has the `medianGQ` and `medianDP` included in the VCF, through the use of `bcftools query`.

## Outcomes

You should now (roughly) know:
- What the requirements are to annotate a VCF
- Know how to annotate a VCF
