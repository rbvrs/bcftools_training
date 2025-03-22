# Annotating your VCF with ClinVar data

For the final exercise, we will annotate our data with [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/) data. ClinVar is a publicly available archive containing information on human genetic variants that are classified for diseases and drug responses, along with supporting evidence. 

When researchers are working with genomic variant data, they will often want to know the clinical impact of a given variant to assess whether it may be causal to the symptoms of a participant (see also [ClinVar's pages on classifications](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) and [the paper on which these are based](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/)). When dealing with very rare variants however, the clinical impact may be unknown. Other databases exist which provide more information on whether a variant could result in a "loss of function" protein, or whether it is predicted to result into a splicing variant causing changes to the resulting transcript.

We already prepared the ClinVar data into a smaller format. The steps we took to generate this data can be found in `/generate_source_material/prepare_clinvar_data.md`.

In this exercise, we will annotate our data with the features listed in the table below (this is just a subset of the data, but these will be sufficient for our purpose). Following the annotation we will be looking at variants of clinical significance and see whether our cohort contains any participants who may carry this variant.

| INFO field | Description |
| ----------- | ----------- |
| CLNDISDB | Tag-value pairs of disease database name and identifier submitted for germline classifications, e.g. OMIM codes and HPO terms indicating human clinical terms |
| CLNDN | ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB |
| CLNHGVS | Variant location and effect. Top-level (primary assembly, alt, or patch) HGVS expression. |
| CLNSIG | Clinical significance. Aggregate germline classification for this single variant; multiple values are separated by a vertical bar |
| GENEINFO | Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|) |


## Exercise requirements

Along with knowledge of some basic linux commands, the following is needed for this exercise.

Required tools:
```
bcftools --version
```

Required (ready made) files:
```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.annotated.filt.vcf.gz"
CLINVAR_DATA="resources/clinvar_20250209.prefixadj.vcf.gz"
```

If you made your own annotation and filter file, adjust the file paths accordingly.

## Exercise

1. Preparation of annotation:
	- The ClinVar data comes with a lot of different INFO fields, a subset was provided in the description above. How many different INFO fields does it contain?
	- Looking back at the previous exercises, what are the three things you will need to annotate your VCF with ClinVar data? 
2. Based on the previous exercises, and the previous question, annotate your VCF with the ClinVar data with the following INFO fields: CLNDISDB, CLNDN, CLNHGVS, CLNSIG, GENEINFO.
3. After annotation, confirm that you have successfully annotated the data:
	- Do you observe the new lines in the header?
	- Using a simply `bcftools query`, can you extract all the possible values for `INFO/CLNSIG`? If so, which one do you observe (optional: and how often do they occur)?
4. After having confirmed that the annotation was performed correctly, we will now look at high quality variants only.
	- Filtering the data for `PASS` variants only, how many variants are classified as `Likely_pathogenic`, and how many as `Uncertain_significance`?
	- How many participants are not homozygote reference for `Likely_pathogenic` variants?
	- How many participants are heterozygote for variants classified as `Uncertain_significance`, and across how many unique variant positions is this?




