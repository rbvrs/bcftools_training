# Implementation of custom filters

The answers are based on the files provided within the repo. If you have files customised yourself, the answers may very depending on the question. In some cases, there are more ways to obtain the answer, the answers provided are one way to do this. Sometimes the answers outlined below are not directly an answer to the questions in the exercise, but you will be able to manually count or infer the answers from them.

**Task: Check previous filters in our VCF**

```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.annotated.vcf.gz"

bcftools view -h $INPUT_FILE | grep "##FILTER"
```

**Task: Introduce filters into our VCF**


```
INPUT_FILE="data/1kgp.chr21.region.sample.subset.annotated.vcf.gz"
ANNOTATION_AND_FILTER_OUTPUT="data/1kgp.chr21.region.sample.subset.annotated.filt.vcf.gz"

bcftools filter --mode + -s low_medianGQ -e 'medianGQ<20' ${INPUT_FILE} -Oz -o tmp/tmp.vcf.gz
bcftools filter --mode + -s low_medianDP -e 'medianDP<10' tmp/tmp.vcf.gz -Oz -o ${ANNOTATION_AND_FILTER_OUTPUT}
bcftools index -t ${ANNOTATION_AND_FILTER_OUTPUT}


bcftools view -H -i 'FILTER="low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
bcftools view -H -i 'FILTER="low_medianDP"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
bcftools view -H -i 'FILTER="low_medianDP;low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l
```

**Task: Confirm functionality**

Checking the automatic integration in the header.
```
bcftools view -h ${ANNOTATION_AND_FILTER_OUTPUT} | grep "^##FILTER="
```

Counting the number of variants with `low_medianGQ` amongst other filters, and uniquely low medianGQ. These should result in 1956 and 772 variants, respectively.
```
bcftools view -H -i 'FILTER~"low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
bcftools view -H -i 'FILTER="low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
```

Counting the number of variants with `low_medianDP`, amongst other filters, and uniquely low medianDP. These should result in 205 and 41 variants, respectively.
```
bcftools view -H -i 'FILTER~"low_medianDP"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
bcftools view -H -i 'FILTER="low_medianDP"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
```

The following variants fail for both failures, amongst other filters, and uniquely this combination. These should result in 159 and 61 variants, respectively.
```
bcftools view -H -i 'FILTER~"low_medianDP;low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
bcftools view -H -i 'FILTER="low_medianDP;low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | wc -l 
```

We can also have a closer look at the variants itself (using uniquely failing for low medianGQ as an example). We will quickly see that these concern multi-allelic variants or MNPs, or are generally seemingly complex variants or in low-complexity regions (poly-T).
```
bcftools view -i 'FILTER="low_medianGQ"' ${ANNOTATION_AND_FILTER_OUTPUT} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/medianGQ\n' | head

bcftools view -i 'FILTER="low_medianDP"' ${ANNOTATION_AND_FILTER_OUTPUT} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/medianDP\n' | head
```
