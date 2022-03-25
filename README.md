# vcf_nirvana_parser
**Parse Nirvana-annotated Somatic SNV/indel vcfs into human-readable table format**
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

| Column  | Explanation  | Link  |   
|---|---|---|
| var_id   | ID of variant (based on position and ref/alt alleles |   |   
| chr  | chromsome of variant  |   |   
| position  | position of variant  |   |   
| ref_allele  | Allele in reference genome  |   | 
| alt_allele  | Alternative allele  |   |   
| variantType  | Type of variant (SNV, insertion, deletion etc)  |   |  
| genotype  | Genotype of site (always 0/1 or 0|1 for somatic)  |   |   
| AF | Alternative Allele Frequency  |   |   
| total_depth  | Approximate read depth at base (reads with MQ=255 or with bad mates are filtered) |   |   
| alt_allele_depth  | Approximate read depth with alternative allele  |   | 
| ref_allele_depth |  Approximate read depth with reference allele  |   |   
| hgnc  |  HUGO Gene Nomenclature ID of gene |   |   
| Ensembl_transcript | Ensembl transcript ID (if variant overlap one such)  | https://www.ensembl.org/index.html  |   
| RefSeq_transcript  | RefSeq transcript ID (if variant overlap one such) | https://www.ncbi.nlm.nih.gov/refseq/  |   
| bioType  | BioType of transcript (protein_coding,lincRNA etc) based on Ensembl/RefSeq annotation  |   |   
| globalMinorAllele  |  dbSNP is reporting the minor allele frequency for each rs included in  a default global population. Since this is being provided to distinguish common polymorphism from rare variants, the MAF is actually the second most frequent allele value. In other words, if there are 3 alleles, with frequencies of 0.50, 0.49, and 0.01, the MAF will be reported as 0.49. The current default global population is 1000Genome phase 3 genotype data from 2500 worldwide individuals, released in the May 2013 dataset.   | https://www.ncbi.nlm.nih.gov/projects/SNP/docs/rs_attributes.html  |   
| globalMinorAlleleFreq  | As above - this is the frequency of GlobalMinorAllele  |   |   
| filterVCdragen  | Filters by Dragen - see explanation in link  | https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/PostSomaticFilters.htm  |   
| SomaticQuality  |  Quality Score of somatic variant call: "You can use somatic quality (SQ) as the primary metric to describe the confidence with which the caller made a somatic call. SQ is a Phred-scaled posterior probability reported as a format field for the tumor sample (exception: for homozygous reference calls in gVCF mode it is instead a likelihood ratio, analogous to hom-ref GQ as described in the germline section). Variants with SQ score below the SQ filter threshold are filtered out using the weak_evidence tag." | https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/SomaticMode.htm  |   
| clinVar_significance  | Clinical Significance of variant (ClinVar)  | ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/intro/, ClinVar Significance https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/  |   
| clinVar_phenotype  | Clinical Phenotypes of variant: ClinVar aggregates the names of medical conditions with a genetic basis from such sources as SNOMED CT, GeneReviews, Genetic Home Reference, Office of Rare Diseases, MeSH, and OMIM®. ClinVar also aggregates descriptions of associated traits from Human Phenotype Ontology (HPO), OMIM, and other sources.   | ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/intro/  |   
| regulatory_region  | Ensembl Regulatory Regions at site  | https://www.ensembl.org/info/genome/funcgen/regulatory_features.html  |   
| revel  | REVEL is an ensemble method for predicting the pathogenicity of missense variants based on a combination of scores from 13 individual tools: MutPred, FATHMM v2.3, VEST 3.0, PolyPhen-2, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP++, SiPhy, phyloP, and phastCons.  |   |   
| dbsnp | dbSNP NCBI annotation of variant  | https://www.ncbi.nlm.nih.gov/snp/ , https://illumina.github.io/NirvanaDocumentation/data-sources/revel |  
| phylopScore  | PhyloP (phylogenetic p-values) conservation scores are obtained from the [PHAST package] (http://compgen.bscb.cornell.edu/phast/) for multiple alignments of vertebrate genomes to the human genome. For GRCh38, the multiple alignments are against 19 mammals and for GRCh37, it is against 45 vertebrate genomes.  | https://illumina.github.io/NirvanaDocumentation/data-sources/phylop	  |   
| primateAI  |  Primate AI is a deep residual neural network for classifying the pathogenicity of missense mutations. The method is described in the publication: | Sundaram, L., Gao, H., Padigepati, S.R. et al. Predicting the clinical impact of human mutation with deep neural networks. Nat Genet 50, 1161–1170 (2018). https://doi.org/10.1038/s41588-018-0167-z  |   
| transcript_consequence  | Ensembl Variation - Calculated variant consequences: "For each variant that is mapped to the reference genome, we identify all overlapping Ensembl transcripts. We then use a rule-based approach to predict the effects that each allele of the variant may have on each transcript. "  | https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html  |   
| aminoAcids  | AminoAcids affected. If missense mutation, it shows RefAA/AltAA. If synonymous mutation, shows the RefAA only  |   |   
| hgvsc  | Nomenclature for variant (cDNA)  |  https://varnomen.hgvs.org/ |   
| hgvsp  | Nomenclature for variant (Protein)  | https://varnomen.hgvs.org/  |   
| polyPhenScore |  PolyPhen2 redicts the effect of an amino acid substitution on the structure and function of a protein using sequence homology, Pfam annotations, 3D structures from PDB where available, and a number of other databases and tools (including DSSP, ncoils etc.). As with SIFT, for each amino acid substitution where we have been able to calculate a prediction, we provide both a qualitative prediction (one of 'probably damaging', 'possibly damaging', 'benign' or 'unknown') and a score. The PolyPhen score represents the probability that a substitution is damaging, so values nearer one are more confidently predicted to be deleterious (note that this the opposite to SIFT). The qualitative prediction is based on the False Positive Rate of the classifier model used to make the predictions.   | http://genetics.bwh.harvard.edu/pph2/dokuwiki/start  | 
| polyPhenPred  | PolyPhen2 (above) predicted effect of mutation on protein function/structure  |   |   
| siftPred  |  SIFT predicts whether an amino acid substitution is likely to affect protein function based on sequence homology and the physico-chemical similarity between the alternate amino acids. The data we provide for each amino acid substitution is a score and a qualitative prediction (either 'tolerated' or 'deleterious'). The score is the normalized probability that the amino acid change is tolerated so scores nearer zero are more likely to be deleterious. The qualitative prediction is derived from this score such that substitutions with a score < 0.05 are called 'deleterious' and all others are called 'tolerated'.  | https://sift.bii.a-star.edu.sg/  |   
| siftScore  |  The score is the normalized probability that the amino acid change is tolerated so scores nearer zero are more likely to be deleterious. |  https://sift.bii.a-star.edu.sg/  |  
| codons  | Altered codon (in case of protein-coding mutations)  |   | 
| exons  | Affected exon (in case of protein-coding mutations)  |   | 
|  introns | Affected intron (in case of intronic mutations) |   |  


																																		<img width="3935" alt="image" src="https://user-images.githubusercontent.com/33324443/160140912-b558f4be-41da-4579-950a-263d6f153d48.png">
