# Exercise cheat sheet

## Table of contents

1. [Introduction](#introduction)
2. [Navigation](#navigation)
3. [Query](#query)
	1. [Querying specific genotypes](#queryA)
4. [Annotation](#annotation)
	1. [Making header entries](#annotationA)
	2. [Preparation of your annotation file](#annotationB)
	3. [Annotating with the file of interest](#annotationC)
5. [Filtering](#filter)

## Introduction <a id="introduction"></a>

`bcftools` is a widely applicable tool to work with variant data, and thus also comes with a [manual](https://samtools.github.io/bcftools/bcftools.html). However, this manual has made no person happy nor has it really clarified anything unless you have already built up some level of perseverance and experience with VCFs and `bcftools` itself anyway. Only then it will start to make more sense, and at that stage, it will become a very powerful tool. The tool is also kept up to date a lot and its community and their developers on [github](https://github.com/samtools/bcftools) are really nice too. 

Once you start to master `bcftools`, it's a very useful skill to have if you intend to work with variant data in your career. We have reported some bugs, or unknown features, to them a few times which they fixed rather quickly. As long as you provide example data for them to test and reproduce the bugs.

Anyway, below we will outline some commands, or parts of it, that you may need during the exercise itself.

## Navigation <a id="navigation"></a>

The below descriptions may not be completely accurate or inclusive, as they may have been shortened in order for them to be applicable to the exercises without causing further complication that the manual itself may bring.

| What? | Command | Note |
| ----------- | ----------- | ----------- |
| View a VCF | `bcftools view <FILE>` | You will also use this to specifically extract parts of the VCF. |
| View a VCF *without* its header. | `bcftools view -H <FILE>` | |
| View a VCF's header | `bcftools view -h <FILE>` | Note, you can pipe the commands into `grep "##contig"` or `grep ##INFO` to capture specific elements. |
| View variant data of a specific sample *without* its header. | `bcftools view -H -s <SAMPLE_NAME> <FILE>` | You can use `-H` or `-h` here, so this is just an example. |
| View variant data of specific samples using an input list *without* its header. | `bcftools view -H -S <SAMPLE_LIST> <FILE>` | You can use `-H` or `-h` here, so this is just an example. |
| View variant data of a specific region *without* its header. | `bcftools view -H -r <LOCATION> <FILE>` | If providing a region, it will follow this notation `CHR:START-END` or `CHR:START`. You can use `-H` or `-h` here, so this is just an example. |
| Output a new VCF | `bcftools view <FILE> -Oz -o <OUTPUT_FILE>` | When generating a new VCF, always aim to maintain the header where possible, hence we do not specify `-H` here. The `-Oz` indicates that the output type needs to be gzipped and the `-o` will be the output file name. |


## Query <a id="query"></a>

The below will outline some of the uses through `bcftools query`, which is another commonly used `bcftools` function. The first table provides an overall outline, where the second table below provides context on the actual query string.

| What? | Command |
| ----------- | ----------- |
| List the samples in a VCF | `bcftools query -l <FILE>` |
| List specific variant info from a VCF | `bcftools query -f '<QUERY_STRING>' <FILE>` |
| List specific variant info of a sample in a VCF | `bcftools query -s <SAMPLE_NAME> -f '<QUERY_STRING>' <FILE>` |
| List specific variant info of a sample in a VCF but only when they are of a specific type (ie. non-ref only, see more options further down) | `bcftools query -s <SAMPLE_NAME> -f '<QUERY_STRING>' -i 'GT="alt"' <FILE>` |

&nbsp;

As mentioned, the `<QUERY_STRING>` above can be quite complex, so we outline a few below that may be useful.

| What? | `<QUERY_STRING>` |
| ----------- | ----------- |
| Extract the genotype `FORMAT` field (`GT`) for each sample | `-f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]'` |
| Extract variant level info, such as `FILTER` | `-f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n'` |
| Extract variant level info, such as an `INFO` field like `AF` (allele frequency), `AC` (allele count), and `AN` (allele number). Note that `INFO` fields always require the following formatting: `%INFO/<name>` | `-f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\n'` |

Note that in the query string's the `[` and `]` are sometimes covering the entire string. This can be indicative that it is *per* sample, and as such you require to place both `\n` within the bracket, as well as outside. Sometimes it may be required to query information that is identical across all samples, so it is variant specific such as the `FILTER` and `INFO` fields. In such cases, you may not want to use those brackets. We recommend playing around with this to see the full effects of it.

### Querying specific genotypes <a id="queryA"></a>

In the first query table, we provided an example where we narrowed down our query to only include variants, leaving out homozygote reference alleles. Below is a table to provide you with other options as well. This table can also be found on the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html), so we simply highlight it here for convenience.

For clarity, if you want to filter your data on genotype (`GT`), you can simply add the following to your bcftools command: `-i 'GT="alt"'` (or any of the other `GT` filters listed below). Note the usage of the various quotation mark types here.


| Filter | Sample genotype |
| ----------- | ----------- |
| `GT="ref"` | reference (haploid or diploid) |
| `GT="alt"` | alternate (hom or het, haploid or diploid) |
| `GT="mis"` | missing genotype |
| `GT="hom"` | homozygous |
| `GT="het"` | heterozygous |
| `GT="hap"` | haploid |
| `GT="RR"` | ref-ref ho |
| `GT="AA"` |  alt-alt hom |
| `GT="RA" or GT="AR"` | ref-alt het |
| `GT="Aa" or GT="aA"` | alt-alt het |
| `GT="R"` | haploid ref |
| `GT="A"` | haploid alt (case-insensitive) |




## Annotation <a id="annotation"></a>

There are two stages to annotating VCFs through `bcftools`. The first stage is generating entries, or lines of text, to add to header, and the second stage is the annotation itself. The example below will focus on INFO field annotation, but `FORMAT` and `FILTER` can be added as well. See sections [1.4.2 to 1.4.4 in the VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) spec for more information.

### Making header entries <a id="annotationA"></a>

Headers will always have a specific format, starting with the `ID`, then `Number` for the number of values within the field, `Type` for the value type, and the general `Description` which is often just freetext. This can be expanded with `Source` and `Version` but these are optional, so we will leave them out for now.

```
touch add_these_lines_to_header.txt
echo '##INFO=<ID=my_annotation_field,Number=1,Type=Float,Description="This value gives context on how much I like this variant.">' >> add_these_lines_to_header.txt
```

### Preparation of your annotation file <a id="annotationB"></a>

If you want to annotate with a `.tsv` file, we recommend gzipping it through `bgzip`. Also make sure that the column names of your file are `CHROM,POS,REF,ALT,my_annotation_field` before running the annotation. After gzipping, index your file through `tabix`.

```
cat my_annotation_file.tsv | bgzip -c > my_annotation_file.tsv.gz
tabix -s1 -b2 -e2 my_annotation_file.tsv.gz
```

### Annotating with the file of interest <a id="annotationC"></a>

Following the steps above, we can now run the annotation as followed. If you want to add more fields to the annotations, you can do this too by just adding more `,+my_second_annotation_field`, provided that you took the necessary steps above too.

```
HEADER_TO_ADD="add_these_lines_to_header.txt"
FILE_WITH_MY_ANNOTATION="my_annotation_file.tsv.gz"
TARGET_FILE_TO_ANNOTATE="target.vcf.gz"
OUTPUT="target.annotated.vcf.gz"

bcftools annotate -h ${HEADER_TO_ADD} \
	-a ${FILE_WITH_MY_ANNOTATION}.gz \
	-c CHROM,POS,REF,ALT,+my_annotation_field ${TARGET_FILE_TO_ANNOTATE} \
	-Oz -o ${OUTPUT}
bcftools index -t ${OUTPUT}
```

## Filtering <a id="filter"></a>

`bcftools` also allows you to build filters. If variants fall within conditions of multiple filters, they will have all of these listed as a `;`-separated `string`, in the `FILTER` field. Anything else will generally be marked as `PASS`.

We should note that when generating a new filter, these unfortunately need to be added one by one.

| What? | Command | Note |
| ----------- | ----------- | ----------- |
| Implement a filter (`-e`) | `bcftools filter --mode + -s my_filter_name -e 'minDP<10' <FILE> -Oz -o <OUTPUT_FILE>` | `--mode +` will indicate that any previous filters are taken along too. bcftools will also generate a line in the header automatically with the `minDP<10` in this case being placed as `Description`. `-e` is exclusionary indicating that all PASS variants will have at least `DP` of 10. |
| Implement a filter (`-i`) | `bcftools filter --mode + -s my_filter_name -i 'minDP<10' <FILE> -Oz -o <OUTPUT_FILE>` | `--mode +` will indicate that any previous filters are taken along too. bcftools will also generate a line in the header automatically with the `minDP<10` in this case being placed as `Description`.  `-i` is inclusive indicating that all PASS variants will have at most `DP` of 10 (generally would not be a filter, so this is illustrative). |
| Viewing all variants of a specific filter | `bcftools view -i 'FILTER="my_filter_name"' <FILE>` |  See also the table below for other expressions to aid filtering. |


After you have generated a custom filter, you can use the expressions below to better filter your data according to your needs. Same as earlier, this (edited) table can also be found on the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html), so we simply highlight it here for convenience.


| Expression | Outcome |
| ----------- | ----------- | 
| FILTER="PASS" | filters for PASS variants only |
| FILTER="." | filters for variants without any flags (often if PASS isn't present) |
| FILTER="A" | exact match, for example "A;B" does not pass |
| FILTER="A;B" | exact match, "A;B" and "B;A" pass, everything else fails |
| FILTER!="A" | exact match, for example "A;B" does pass |
| FILTER~"A" | subset match, for example both "A" and "A;B" pass |
| FILTER~"A;B" | subset match, pass only if both "A" and "B" are present |
| FILTER!~"A" | complement match, for example both "A" and "A;B" fail |
| FILTER!~"A;B" | complement match, fail if both "A" and "B" are present | 


