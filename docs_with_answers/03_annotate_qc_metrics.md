# Annotation of QC metrics

The answers are based on the files provided within the repo. If you have files customised yourself, the answers may very depending on the question. In some cases, there are more ways to obtain the answer, the answers provided are one way to do this. Sometimes the answers outlined below are not directly an answer to the questions in the exercise, but you will be able to manually count or infer the answers from them.

**Task: Generate new header lines**

Make a header file containing medianGQ and medianDP information. These are `INFO` fields because they are variant-level data.

```
ANNOTATION_HEADER="aux/add_qc_metrics_info_lines_to_header.txt"

touch ${ANNOTATION_HEADER}

echo '##INFO=<ID=medianGQ,Number=1,Type=Float,Description="Median genotype quality across samples">' > ${ANNOTATION_HEADER}
echo '##INFO=<ID=medianDP,Number=1,Type=Float,Description="Median depth across samples">' >> ${ANNOTATION_HEADER}
```

**Task: Gzip and index the annotation file**

Annotation is significantly made faster and easier for `bcftools` when the file we work with is indexed. We therefore also ensure that the file itself is first gzipped and indexed with `bgzip` and `tabix` respectively.

```
ANNOTATION_FILE="aux/gq_dp_annotation_file.txt"

cat ${ANNOTATION_FILE} | bgzip -c > ${ANNOTATION_FILE}.gz
tabix -s1 -b2 -e2 ${ANNOTATION_FILE}.gz
```

**Task: Annotate the VCF with QC metrics**

We are first setting all the variables and datasets that we require. Once set we then run the annotation itself.

```
ANNOTATION_HEADER="aux/add_qc_metrics_info_lines_to_header.txt"
ANNOTATION_FILE_ZIPPED="aux/gq_dp_annotation_file.txt.gz"
TARGET_FILE="data/1kgp.chr21.region.sample.subset.vcf.gz"

ANNOTATED_VCF="data/1kgp.chr21.region.sample.subset.annotated.vcf.gz"

bcftools annotate -h ${ANNOTATION_HEADER} \
	-a ${ANNOTATION_FILE_ZIPPED} \
	-c CHROM,POS,REF,ALT,+medianGQ,+medianDP ${TARGET_FILE} \
	-Oz -o ${ANNOTATED_VCF}
bcftools index -t ${ANNOTATED_VCF}
```

**Task: Confirming that medianGQ and medianDP have been integrated**

After annotation we can quickly check whether the annotation worked. First we provide a line to check the contents of the header. And afterwards, we run a `bcftools query` to specifically extract the `medianGQ` and `medianDP` values again. When the annotation doesn't work, the `bcftools query` should result in an error stating that one of the `INFO` fields is missing.

```
ANNOTATED_VCF="data/1kgp.chr21.region.sample.subset.annotated.vcf.gz"

bcftools view -h ${ANNOTATED_VCF} | grep -E "medianDP|medianGQ"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/medianGQ\t%INFO/medianDP\n' ${ANNOTATED_VCF} | head
```
