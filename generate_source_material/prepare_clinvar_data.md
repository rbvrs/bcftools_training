# ClinVar

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/) is a publicly available archive containing information on human genetic variants that are classified for diseases and drug responses, along with supporting evidence.

We will be using this data as an example on how variant level data can be annotated and filtered for potential pathogenic variants of interest within a cohort of individuals. As an example, this data will contain a column (`CLINSIG`) indicating clinical significance, on which more information can be found [here](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/).

For this exercise we will focus on a specific region within chromosome 21 to maintain low computational requirements, however it can be expanded to any other region or chromosome of interest.

## Download ClinVar Data

ClinVar provides its data in a clean VCF format, and gets updated every single [week](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/). It is therefore advisable to use the latest version to work with where possible. The main VCF (`clinvar_yyyymmdd.vcf.gz`) is around 140-150 Mb.

Please find the ftp link for the latest version [here](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/). We will be working with variant data based on the `GRCh38` refrence. The file we will be using is from the 9th of February, 2025: `clinvar_20250209.vcf.gz`. If you wish to download a later version, it should not significantly affect the downstream material, though the numbers *may* be slightly off.

Please use the below code snippet to download your data into the `resources` folder.
```
URL_PATH="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/"
FILE="clinvar_20250209.vcf.gz"

echo "Downloading from: ${URL_PATH}${FILE}"
curl -O ${URL_PATH}${FILE}  --output-dir "resources"
curl -O ${URL_PATH}${FILE}.tbi  --output-dir "resources"
```
You should now have the ClinVar data at `resources/clinvar_20250209.vcf.gz`.

## Adjust chromosome prefix (required!)

[VCFs](https://samtools.github.io/hts-specs/VCFv4.3.pdf) will always have a specific tab-delimited notation in its body. Often this is `CHROM-POS-REF-ALT` as separate columns indicating the location of a variant. As mentioned, we will be working with variants aligned to the GRCh38 reference genome.

You may learn that, depending on the reference that was used, chromosome prefixes may be absent or present. Generally, for GRCh38 chromosomes will have a `chr` prefix like `chr1`, `chr17`, or `chrX`. However, in other cases there may only be `1`, `17`, or `X`. ClinVar will not have the `chr` prefix, but the data that we will be working with, will have this prefix and therefore we need to adjust the data. Note, that just changing this prefix will not make data compatible from GRCh37 to GRCh38 or vice versa given positions between these two references are different. 

To add these prefixes, we will use a mapping file (`aux/grch38_chr_mappings.txt`) and `bcftools` which should be part of the conda environment.


```
CLINVAR_INPUT="resources/clinvar_20250209.vcf.gz"
MAPPING_FILE="aux/grch38_chr_mappings.txt"

OUTPUT_FILE="resources/clinvar_20250209.prefixadj.vcf.gz"

bcftools annotate --rename-chrs ${MAPPING_FILE} ${CLINVAR_INPUT} -Oz -o ${OUTPUT_FILE}
bcftools index -t ${OUTPUT_FILE}
```


## Confirm `chr` prefix adjustment

Run the following command to check the changes in the body (note `-H`) of the VCF:

`bcftools view -H ${OUTPUT_FILE} | head -n 2`.

This should result in the following:
```
chr1    66926   3385321 AG      A       .       .       ALLELEID=3544463;CLNDISDB=Human_Phenotype_Ontology:HP:0000547,MONDO:MONDO:0019200,MeSH:D012174,MedGen:C0035334,OMIM:268000,OMIM:PS268000,Orphanet:791;CLNDN=Retinitis_pigmentosa;CLNHGVS=NC_000001.11:g.66927del;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGSCV=SCV005419006;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=OR4F5:79501;MC=SO:0001627|intron_variant;ORIGIN=0
chr1    69134   2205837 A       G       .       .       ALLELEID=2193183;CLNDISDB=MedGen:CN169374;CLNDN=not_specified;CLNHGVS=NC_000001.11:g.69134A>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGSCV=SCV003526545;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=ClinGen:CA502008;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1```
```

Second, we will confirm the changes in the header (note `-h`) of the VCF. Run the following command: 

`bcftools view -h ${OUTPUT_FILE} | grep "##contig"`

This should result in the following:
```
##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>
##contig=<ID=chr6>
##contig=<ID=chr7>
##contig=<ID=chr8>
##contig=<ID=chr9>
##contig=<ID=chr10>
##contig=<ID=chr11>
##contig=<ID=chr12>
##contig=<ID=chr13>
##contig=<ID=chr14>
##contig=<ID=chr15>
##contig=<ID=chr16>
##contig=<ID=chr17>
##contig=<ID=chr18>
##contig=<ID=chr19>
##contig=<ID=chr20>
##contig=<ID=chr21>
##contig=<ID=chr22>
##contig=<ID=chrX>
##contig=<ID=chrY>
##contig=<ID=chrM>
##contig=<ID=NT_113889.1>
##contig=<ID=NT_187633.1>
##contig=<ID=NT_187661.1>
##contig=<ID=NT_187693.1>
##contig=<ID=NW_009646201.1>
```

You can run the same commands for the original file too to get a better idea of the changes.

You should now be ready to use the adjusted ClinVar data for the rest of the exercises.