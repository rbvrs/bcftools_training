# README

## Background

This repo contains exercises focussed around the use of `bcftools` and various sub functions like `view`, `query`, and `annotate`. As an annotation example, we use ClinVar to annotate our data.

The exercises are divided in order to focus on specific functionality. A parallel directory is provided with answers and explanations. However, we also provide a cheat-sheet which we encourage to be use alongside the regular exercises as it covers all the functionality that is required within the exercises.

## Required skillset

Exposure to command line interface will be needed to successfully complete these exercises. Beyond that, the majority of exercises will be `bcftools` commands which will be outlined step by step throughout the material along with a 'cheat sheet'. 
No starting knowledge of `bcftools` is required as that is the purpose of the material. 

Some knowledge on linux commands (e.g. `grep`, `sort`, `cut`, `uniq`, `wc -l`) will be needed too. The requirement of `awk` is avoided to keep the requirements accessible, but this can be used to complete the exercises too.

The exercises are run from a conda environment, so basic knowledge on initiating a conda environment will be needed too.

`Python` or `R` data wrangling skills may be useful, but not a hard requirement, for exercise `02`. The intermediate files are provided in case this skillset is not yet developed.

## Environment

The tools required for these training exercises can be run through a conda environment. The `.yaml` file for this environment is located in `/conda-env/training-env.yaml`.

### Required tools within environment

The following tools should be part of your environment to run these exercises successfully:

```
bcftools --version
bgzip --version
tabix --version
```

### Download intermediate files

After cloning the repo to your local machine or compute system, you can download the following intermediate files required throughout the exercises from [Zenodo](https://zenodo.org/records/15068836). The total size of this data altogether is ~3.5 Gb. The largest file is a `.tsv` which will be located in the `aux/` folder, and sits at 2 Gb. The remainder are VCFs ~260 Mb max.

You can download all the required files through the following `wget` commands.

```
wget -P data/ -i aux/download_data.txt
wget -P tmp/ -i aux/download_tmp.txt
wget -P resources/ -i aux/download_resources.txt
wget -P aux/ -i aux/download_aux.txt
```

## Source material and inspiration

### Source data

The exercises made available in this repo are based on [1KGP](https://www.internationalgenome.org/1000-genomes-summary) data. Please cite these accordingly if these are further used.
For clinical annotation data, we work with [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/). Again here, please do cite them accordingly.
Please find additional documentation on the source material in [generating source data](generate_source_material/).

### Material inspiration and other resources

The majority of the material presented here is written and designed by myself, but I have used the following resources as conceptual inspiration.

- [Glittr](https://glittr.org/?per_page=25&sort_by=stargazers&sort_direction=desc) (Git repositories with educational materials for the computational life sciences). I want to particularly highlight this as a good resource to find further training material such as the Harvard Chan Bioinformatics Core lecture material (below). Glittr is built and maintained by the [Swiss Institute of Bioinformatics](https://www.sib.swiss/). I highly recommend any lecturer or course organiser to browse through this resource.
- [Harvard Chan Bioinformatics Core lecture material on variant analysis](https://github.com/hbctraining/Intro-to-variant-analysis/tree/main/lectures).
- [VCF format PDF](https://samtools.github.io/hts-specs/VCFv4.3.pdf), as well as the following [GATK resource](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format) on VCF.
- [BCFtools manual](https://samtools.github.io/bcftools/bcftools.html). Great manual, but complex so it has quite a high entry level.

### Disclaimer

This training material is not formally endorsed, supported, or maintained by [Genomics England Ltd.](https://www.genomicsengland.co.uk/) and was made in my own name during my employment with them. No data from Genomics England Ltd. is used during these exercise.
While I do not aim to expand on the material provided here, feel free to raise any issues/questions/feedback directly on the repo. 
I will try to attend to them where appropriate. 

