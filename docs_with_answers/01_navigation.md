# Navigation answers

The answers are based on the files provided within the repo. If you have files customised yourself, the answers may very depending on the question. In some cases, there are more ways to obtain the answer, the answers provided are one way to do this. Sometimes the answers outlined below are not directly an answer to the questions in the exercise, but you will be able to manually count or infer the answers from them.

1. 118,302 variants
	- `bcftools -H view ${INPUT_FILE} | wc -l`
2. 12x `FORMAT` fields, 95x `INFO` fields
	- `bcftools -H view ${INPUT_FILE} | grep "##FORMAT" | wc -l`
	- `bcftools -H view ${INPUT_FILE} | grep "##INFO" | wc -l`
3. 3x for `HG00330`, and we observe MNPs and multi-allelic variants.
	- `bcftools query -s "HG00330,HG02146" -r "chr21:25775000-25776000" -i 'GT="alt"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' ${INPUT_FILE}`
	- `bcftools query -s "HG00330,HG02146" -r "chr21:25775000-25776000" -i 'GT="alt"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' ${INPUT_FILE} | grep "HG00330" | wc -l`
4. `AF`, `AF_EUR`, `AF_AFR`
	- `bcftools query -r "chr21:27244663" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_EUR\t%INFO/AF_AFR\n' ${INPUT_FILE}`: 0.239382, 0.0315956, 0.681411
	- `bcftools query -r "chr21:25763936" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_EUR\t%INFO/AF_AFR\n' ${INPUT_FILE}`: 0.413043, 0, 0.415
