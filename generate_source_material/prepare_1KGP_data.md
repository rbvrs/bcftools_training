# 1000 Genomes Data

The [1000 Genomes Project](https://www.internationalgenome.org/1000-genomes-summary) (or 1KGP) aimed to discover common variants of at least 1% (`AF>0.01`) across populations. Originally this started with 1,000 genomes as the costs of sequencing dropped further, but has later been amended with additional genomes. This data serves as a good high quality reference set for population genomics, as well as for testing existing pipelines that researchers may have built. See also their main papers on the dataset itself on [human genetic variation](https://www.nature.com/articles/nature15393), and one specifically for [structural variation](https://www.nature.com/articles/nature15394).

As many cohorts and consortia maintain their datasets behind information governance and data protection regulations, the 1KGP dataset helps us to still learn and educate people on working with large variant based datasets.

## Download 1KGP data

The 1KGP offers many various formats of their [data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/), for example a phased dataset will only contain genotype (`GT`) along with a few INFO fields related to allele counts amongst frequency others. However, for the purpose of our exercise, we will work with the recalibrated variants (a slightly more raw version) dataset. These VCFs will contain genotype, allele depth, read depth, and genotype quality as part of their FORMAT fields.

The data can be quite large, so we will initially be focussing on chromosome 21. Even still, this file is 24 Gb so we will be downscaling it to a specific region. Because of its size, we will not place this file on our repo, but below we provide the steps towards downloading it.

The file with chromosome 21 data can be downloaded as followed:

```
URL_PATH="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/"
FILE="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"

curl -O ${URL_PATH}${FILE}  --output-dir "resources"
curl -O ${URL_PATH}${FILE}.tbi  --output-dir "resources"
```

## Confirm file integrity

You can count the number of samples within the file as a quick way to confirm whether the download was successful. The 1KGP dataset also contains a file with `md5sums` (`20201028_genotyped-manifest.tsv`), but that is outside of the scope of this exercise.

The below code snippet should result in a count of `3202` samples.
```
CHR21_1KGP="resources/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"

bcftools query -l ${CHR21_1KGP} | wc -l
```

## Downsize file

There are various ways of downsizing VCFs, either by region or subsetting them by samples. In this case, we will already provide you with a downsized file based on region. For reference, the original file was 24GB containing 3202 samples. The below numbers may vary on the region or samples, so take it as an approximate.

- By region: 	`1.6 Gb`
- By samples:	`4.3 Gb`
- By region + samples:	`260 Mb`

### By region 

In our exercise we will work with a specific region, which contains the  [APP](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000142192;r=21:25880535-26171128) gene. Variation within APP (amyloid precursor protein) is often implicated in Alzheimer's disease (see [Hampel et al., 2021](https://www.nature.com/articles/s41380-021-01249-0#Sec3), [Sirisi et al., 2024](https://alzres.biomedcentral.com/articles/10.1186/s13195-024-01504-w)).

Furthermore, the germline aggregate of Genomics England ([AggV2](https://re-docs.genomicsengland.co.uk/aggv2/)) contains variant data which have been split into 1,371 chunks. This allows for faster querying of individual regions. For our exercise, we will use an identical chunk size for chromosome 21 in which APP is located.

Splitting a VCF into a specific region can be done with the code snippet below:

```
CHUNK_REGION_WITH_APP="chr21:24481187-27270147"
DOWNSIZED_VCF="data/1kgp.chr21.region.vcf.gz"

bcftools view -r ${CHUNK_REGION_WITH_APP} ${CHR21_1KGP} -Oz -o ${DOWNSIZED_VCF}
bcftools index -t ${DOWNSIZED_VCF}
```


### By sample

To subset the data by sample list, we run the following code. We first retrieve a list of all the samples, before randomly selecting 500 samples and placing this in a separate file. Once we have selected 500 random samples, we subset the vcf through the use of `-S`. If we want to select one sample, you can use `-s <SAMPLE_NAME>`, or even multiple samples through comma separation. But once you are working with more samples, it is easier to simply use a separate file.

Generate a list of 500 randomly selected samples:

```
TMP_LIST="tmp/sample_list.txt"
LIST_500="aux/sample_list_random500.txt"

bcftools query -l ${CHR21_1KGP} > ${TMP_LIST}

awk 'BEGIN {srand()} {print rand(), $0}' ${TMP_LIST} | sort -k1,1n | cut -d' ' -f2- | head -n 500 | sort > ${LIST_500}
```

Generate a new VCF containing only those 500 samples:
```
DOWNSIZED_VCF="data/1kgp.chr21.subset.vcf.gz"

bcftools view -S ${LIST_500} ${CHR21_1KGP} -Oz -o ${DOWNSIZED_VCF}
bcftools index -t ${DOWNSIZED_VCF}
```

### By both sample and region in one command

We can do the above in one go as well. We will assume that you have generated a randomised list of 500 samples.

```
LIST_500="aux/sample_list_random500.txt"
CHUNK_REGION_WITH_APP="chr21:24481187-27270147"
CHR21_1KGP="resources/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz"

OUTPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"

bcftools view -r ${CHUNK_REGION_WITH_APP} -S ${LIST_500} ${CHR21_1KGP} -Oz -o ${OUTPUT_FILE}
bcftools index -t ${OUTPUT_FILE}
```

This file should be approximately 260 Mb, but we can run some other quick checks to see if everything worked out fine. First, we can check the number of samples in the final file, the following line should result in `500` being returned.

```
OUTPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
bcftools query -l ${OUTPUT_FILE} | wc -l
```

We can further count the number of variant locations that are part of the body of the VCF through the following command (answer: `118032`):

```
OUTPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
bcftools view -H ${OUTPUT_FILE} | wc -l 
```

We previously did not run this command on the large VCF as it may take a long time depending on your compute. We can however check the first 5 lines and last 5 lines of the file through the following lines. Note that this will likely flood your terminal, because of all the sample information for each of the 500 samples. However, you should see positions such as `chr21:24481205` through `chr21:24481304` in the `head` and `chr21:27269763` through `chr21:27270066` in the `tail` section.

```
OUTPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
bcftools view -H ${OUTPUT_FILE} | head
bcftools view -H ${OUTPUT_FILE} | tail
```

If you want to do this without (relatively speaking) flooding the terminal, we can also just select the first sample and run the commands above again to account for it.
```
OUTPUT_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"
EXAMPLE_SAMPLE="HG00096"

bcftools view -H -s ${EXAMPLE_SAMPLE} ${OUTPUT_FILE} | head
bcftools view -H -s ${EXAMPLE_SAMPLE} ${OUTPUT_FILE} | tail
```
We should now be ready to work through the rest of the exercises.


