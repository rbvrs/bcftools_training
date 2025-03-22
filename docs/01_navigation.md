# VCF Navigation

This section has a few questions to help you navigate within a VCF.

## Exercise requirements

Along with knowledge of some basic linux commands, the following is needed for this exercise.

Required tools:
```
bcftools --version
```

Required (ready made) files:
```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
```

## Questions

1. How many variants are there within this file?
2. `FORMAT` fields are sample-specfic, whereas `INFO` fields are variant-specific. How many `FORMAT` fields are provided in the VCF? And `INFO` fields?
3. Obtain all heterozygote genotypes for two samples (i.e. `HG00330` and `HG02146`) in a 5000-bp range  starting from `chr21:25775000`. How many unique variants do you see for sample `HG00330`? And what other type of variants do you observe beyond regular SNPs?
3. What is the allele frequency of variant `chr21:27244663`? And what about the allele frequency's within individuals with European ancestry, or with African ancestry? And what about variant `chr21:25763936`? What does this tell you?


## Outcomes

You should now (roughly) know:
- How to navigate through a VCF
- How and where to find information on FORMAT or INFO fields
- Obtain sample-specific (`FORMAT`) information across a list of variation
- Obtain variant-specific (`INFO`) information for specific variants.