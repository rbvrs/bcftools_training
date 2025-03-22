# Add custom filters to our annotated VCF

This section is directed at adding custom filters to our QC-annotated VCF from the previous exercise. Researchers will often want to work with high quality variants. VCFs will have a specific column included, called the `FILTER` column, containing quality information on variants. If a variant passes all the filters without issue, the `FILTER` column will state `PASS` (a VCF may be provided without PASS filter however). If not listed with a `PASS` filter or if it contains a `.`, these filters will indicate whether that variant has particular undesirable quirks or features. Sometimes a variant will have multiple characteristics flagged in case they fail two or mmore filters. In such cases the `FILTER` will state `odd_characteristic_1;odd_characteristic_5` instead of `PASS`. In this exercise we will implement two filters, one for medianGQ and another one for medianDP.

## Exercise requirements

Along with knowledge of some basic linux commands, the following is needed for this exercise.

Required tools:
```
bcftools --version
```

Required (ready made) files:
```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.annotated.vcf.gz"
```
If you made your own annotation file, adjust the file paths accordingly.

## Exercise

Our data will already have some filters included. We will try to maintain but ignore these.

1. Which `FILTER`s are already included in our VCF before we start? You don't necessarily need to know what they do at this stage.
2. Generate filters for medianGQ and medianDP.
	- The medianGQ filter should be set to FAIL for variants that have a medianGQ<20.
	- The medianDP filter should be set to FAIL for variants that have a medianDP<10.
	- Confirm that the filters have indeed been implemented correctly by looking up a variant that should PASS and one that should FAIL.
3. Are there any variants which fail both filters?
4. What do you notice about the variants that fail for these filters?

## Outcomes

You should now (roughly) know:
- What a FILTER is, and how it can be used
- Know how to FILTER a VCF based on metrics
- Know how to build a custom FILTER
