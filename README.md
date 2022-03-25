# vcf_nirvana_parser
**Parse Nirvana-annotated vcfs into human-readable table format**
* NB! Only tested / developed for SNV Dragen VCFs annotated with Nirvana GRCh37'
* Takes only Canonical transcripts (Ensembl or RefSeq) (Future versions will have option to use all transcripts - canonical or not)

#### Nirvana Annotator
https://illumina.github.io/NirvanaDocumentation/

## Input
* Nirvana snv vcf output (.json.gz) format. All samples must be in the same folder (specified by --nirvpath)
e.g. if running with `--nirvpath ./nirvana_files` you need the following structure
```
- nirvana_files
|_sample1.hard-filtered.annotated.nirvana.txt.json.gz
|_sample2.hard-filtered.annotated.nirvana.txt.json.gz
```


## Output
* One table pr sample snv vcf (`SAMPLE1.snv.annotation.nirvana.table.csv`)
* One "PASS" joint table with all samples (tables above) joined, where each sample is represented by AF: only passed variants (`joint_variant_AF_table.PASS.csv`)
* One "Filter" joint table with all samples (tables above) joined, where each sample is represented by AF: all variants, those filtered by Dragen will have an (f) flag behind AF, and an additional column pr sample with the filter_name (`joint_variant_AF_table.FILTER.csv`)

## How to run
`python vcf_parse_to_table.py --nirvpath /path/to/nirvana-dir`

## Output explanation

| Column  | Explanation  | Link  |   |   |
|---|---|---|---|---|
| var_id   | ID of variant (based on position and ref/alt alleles |   |   |   |
| chr  | chromsome of variant  |   |   |   |
| position  | position of variant  |   |   |   |
| ref_allele  | Allele in reference genome  |   |   |   |
| alt_allele  | Alternative allele  |   |   |   |
| variantType  | Type of variant (SNV, insertion, deletion etc)  |   |   |   |
| genotype  | Genotype of site (always 0/1 or 0|1 for somatic)  |   |   |   |
| AF | Alternative Allele Frequency  |   |   |   |
| total_depth  | Approximate read depth at base (reads with MQ=255 or with bad mates are filtered) |   |   |   |
| alt_allele_depth  | Approximate read depth with alternative allele  |   |   |   |
| ref_allele_depth |  Approximate read depth with reference allele  |   |   |   |
| hgnc  |  HUGO Gene Nomenclature ID of gene |   |   |   |
| Ensembl_transcript | Ensembl transcript ID (if variant overlap one such)  | https://www.ensembl.org/index.html  |   |   |
| RefSeq_transcript  | RefSeq transcript ID (if variant overlap one such) |   |   |   |
| bioType  | BioType of transcript (protein_coding,lincRNA etc) based on Ensembl/RefSeq annotation  |   |   |   |
| globalMinorAllele  |   dbSNP is reporting the minor allele frequency for each rs included in  a default global population. Since this is being provided to distinguish common polymorphism from rare variants, the MAF is actually the second most frequent allele value. In other words, if there are 3 alleles, with frequencies of 0.50, 0.49, and 0.01, the MAF will be reported as 0.49. The current default global population is 1000Genome phase 3 genotype data from 2500 worldwide individuals, released in the May 2013 dataset.   | https://www.ncbi.nlm.nih.gov/projects/SNP/docs/rs_attributes.html  |   |   |
| globalMinorAlleleFreq  | As above - this is the frequency of GlobalMinorAllele  |   |   |   |
| filterVCdragen  | Filters by Dragen - see explanation in link  | https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/PostSomaticFilters.htm  |   |   |
| SomaticQuality  |  Quality Score of variant call |   |   |   |
| clinVar_significance  |   |   |   |   |
| clinVar_phenotype  |   |   |   |   |
| regulatory_region  |   |   |   |   |
| revel  |   |   |   |   |
| dbsnp |   |   |   |   |
| phylopScore  |   |   |   |   |
| primateAI  |   |   |   |   |
| transcript_consequence  |   |   |   |   |
| aminoAcids  |   |   |   |   |
| hgvsc  |   |   |   |   |
| hgvsp  |   |   |   |   |
|  polyPhenScore |   |   |   |   |
| polyPhenPred  |   |   |   |   |
| siftPred  |   |   |   |   |
| siftScore  |   |   |   |   |
| codons  |   |   |   |   |
| exons  |   |   |   |   |
|  introns |   |   |   |   |
|   |   |   |   |   |


																																		<img width="3935" alt="image" src="https://user-images.githubusercontent.com/33324443/160140912-b558f4be-41da-4579-950a-263d6f153d48.png">
