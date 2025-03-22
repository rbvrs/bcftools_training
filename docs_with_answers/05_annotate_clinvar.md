# Annotation of your VCF with ClinVar data

The answers are based on the files provided within the repo. If you have files customised yourself, the answers may very depending on the question. In some cases, there are more ways to obtain the answer, the answers provided are one way to do this. Sometimes the answers outlined below are not directly an answer to the questions in the exercise, but you will be able to manually count or infer the answers from them.

**Task: Annotation preparation of ClinVar data**

Note that we can now annotate our VCF directly with the ClinVar VCF, and that we do not need to extract the individual INFO fields through `bcftools query`.

Using the following command we find how many INFO fields are contained. We specifically look at the header of the ClinVar VCF and `grep` out the specific `##INFO` term to count how many INFO fields there are present.

```
CLINVAR_DATA="resources/clinvar_20250209.prefixadj.vcf.gz"
bcftools view -h ${CLINVAR_DATA} | grep "##INFO" | wc -l
```
Answers:
- The number of INFO fields in the VCF is `39`
- An annotation header with all the descriptions for the INFO fields of interest, an input file (the ClinVar VCF), and the file we want to use to annotate our data with 

**Task: Annotation with ClinVar data**

First, we will produce a new header. We do not need to write custom fields in our case, because we can just take the entries from the source ClinVar VCF. We know how to specifically look at the header of the ClinVar data, and how to extract information from it. In this case, we will use a single `grep` command to target and extract elements of the ClinVar header. By using `grep -E` we signify that multiple terms (separated by `|`) can be extracted. We added a `,` after each term because in the header itself it will indicate the end of the term. This may not always be a robust approach, but will suffice for the purpose of our exercise.

```
CLINVAR_DATA="resources/clinvar_20250209.prefixadj.vcf.gz"
CLINVAR_HEADER_INFO="aux/add_clinvar_annotations_to_header.txt"

bcftools view -h ${CLINVAR_DATA} | grep -E "CLNDISDB,|CLNDN,|CLNHGVS,|CLNSIG,|GENEINFO," > ${CLINVAR_HEADER_INFO}
```

This should now result in a header with 5 lines, one per `INFO` field of interest. We have our header, our ClinVar VCF as input data, and the target VCF of interest.

```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.annotated.filt.vcf.gz"
CLINVAR_DATA="resources/clinvar_20250209.prefixadj.vcf.gz"
CLINVAR_HEADER_INFO="aux/add_clinvar_annotations_to_header.txt"

CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"

bcftools annotate -h ${CLINVAR_HEADER_INFO} -a ${CLINVAR_DATA} \
	-c CHROM,POS,REF,ALT,+CLNDISDB,+CLNDN,+CLNHGVS,+CLNSIG,+GENEINFO ${INPUT_FILE} \
	-Oz -o ${CLINVAR_ANNOTATED_OUTPUT}
bcftools index -t ${CLINVAR_ANNOTATED_OUTPUT}
```

**Task: Confirm and Filtering of ClinVar annotation**

The `${CLINVAR_ANNOTATED_OUTPUT}` should now contain the newly annotated data. This can be confirmed in the header through the line below where we should see the newly added fields at the bottom of the header.

```
CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"
bcftools view -h ${CLINVAR_ANNOTATED_OUTPUT} | grep "##INFO"
```

Alternatively, we can specifically use bcftools query to target one of the `INFO` fields, and extra the unique values within it.

```
CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"
bcftools query -f '%INFO/CLNSIG\n' ${CLINVAR_ANNOTATED_OUTPUT} | sort -u
```

However, to count the number of occurences of each classification, run the following line.

```
CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"
bcftools query -f '%INFO/CLNSIG\n' ${CLINVAR_ANNOTATED_OUTPUT} | sort | uniq -c
```

Finally, we count the number of variants that are PASS to ensure that we lower the chance that low quality variants obscure our assessment. On a quick glance, it looks like only four variants were dropped which were potentially related to clinical findings. The majority (~8,000) of variants that were dropped were not present in the ClinVar data, or did not have a clinical significance.

```
CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"
bcftools query  -i 'FILTER="PASS"' -f '%INFO/CLNSIG\n' ${CLINVAR_ANNOTATED_OUTPUT} | sort | uniq -c
```

To directly calculate the number of likely pathogenic variants while also extracting the genotypes of the participants, we can run the following lines of code. The first `bcftools query` filters all the variants for `PASS` only and variants that have `Likely_pathogenic` annotated within the `CLNSIG` `INFO` field. The result of this is a long list, per sample, of all genotypes for these variants. This will thus also include reference genotypes. These can be filtered out by including `GT!="RR"` in the filters. Running this will likely result in 0 lines, because we do not have any participant in our dataset who is non-Ref for this variant.

```
CLINVAR_ANNOTATED_OUTPUT="data/1kgp.chr21.region.sample.subset.clinvar.annotated.filt.vcf.gz"

bcftools query -i 'FILTER="PASS" && INFO/CLNSIG="Likely_pathogenic"' -f "[%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%SAMPLE\t%GT\n]" ${CLINVAR_ANNOTATED_OUTPUT} | head

bcftools query -i 'FILTER="PASS" && INFO/CLNSIG="Likely_pathogenic" && GT!="RR"' -f "[%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%SAMPLE\t%GT\n]" ${CLINVAR_ANNOTATED_OUTPUT} | head
```

As an added example, we can look at variants of uncertain significance (VUS). This can happen when ["a variant does not fulfill criteria using either of these sets (pathogenic or benign), or the evidence for benign and pathogenic is conflicting, the variant defaults to Uncertain Significance"](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/). Filtering for these type of variants, we do retrieve multiple samples that are at least heterozygote for some of them. 

```
bcftools query -i 'FILTER="PASS" && INFO/CLNSIG="Uncertain_significance" && GT!="RR"' -f "[%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%SAMPLE\t%GT\n]" ${CLINVAR_ANNOTATED_OUTPUT} | head
```
To assess and count whether other genotypes are present in the output we can run the following, where we essentially drop the sample column, and count the number of occurrences of each row. The majority of heterozygote variants are only observed in one sample, along with three other variants that are present in two to four samples.

```
bcftools query -i 'FILTER="PASS" && INFO/CLNSIG="Uncertain_significance" && GT!="RR"' -f "[%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%SAMPLE\t%GT\n]" ${CLINVAR_ANNOTATED_OUTPUT} | cut -d$'\t' -f1-4,7 | uniq -c
```