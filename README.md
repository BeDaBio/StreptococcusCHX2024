# Genomic and Transcriptomic Insights into Chlorhexidine Adaptation in Oral Streptococcus spp.: Implications for Antimicrobial Resistance 
Bernd Daller1, David L. Auer2, Wolfgang Buchalla2, Sibylle Bartsch3, André Gessner1, Nicholas S. Jakubovics4, Ali Al-Ahmad3, Andreas Hiergeist1,#, Fabian Cieplik3,#,* 

 

1 Institute of Clinical Microbiology and Hygiene, University Hospital Regensburg, Regensburg, Germany 

2 Department of Conservative Dentistry and Periodontology, University Hospital Regensburg, Regensburg, Germany 

3 Department of Operative Dentistry and Periodontology, Center for Dental Medicine, Medical Center – University of Freiburg, Medical Faculty, University of Freiburg, Freiburg i. Br., Germany 

4 School of Dental Sciences, Faculty of Medical Sciences, Newcastle University, Newcastle upon Tyne, UK 


------------------------------------------------------

This study investigated CHX adaptation in three Streptococcus spp. clinical isolates and their CHX-adapted counterparts using whole genome sequencing and transcriptomic profiling.

Follow these steps to replicate the RNAseq Analysis:

### 1. Clone repository containing relevant Data:
```bash
Data
├── Antibiotic_Resistance
│   └── hamronization_combined_report.tsv
├── Bakta_annotated_files
│   ├── Strain_73
│   │   ├── Strain_73_wt.fna
│   │   ├── Strain_73_wt.gff3
│   │   └── Strain_73_wt.tsv
│   ├── Strain_78
│   │   ├── Strain_78_wt.fna
│   │   ├── Strain_78_wt.gff3
│   │   └── Strain_78_wt.tsv
│   └── Strain_93
│       ├── Strain_93_wt.fna
│       ├── Strain_93_wt.gff3
│       └── Strain_93_wt.tsv
├── GO_Check
│   ├── go-basic.obo
│   └── gocheck_do_not_annotate.obo
├── GO_Pannzer2
│   ├── DE_strain_73_wt.csv
│   ├── DE_strain_78_wt.csv
│   ├── DE_strain_93_wt.csv
│   ├── GO_strain_73_wt.csv
│   ├── GO_strain_78_wt.csv
│   ├── GO_strain_93_wt.csv
│   ├── strain_73_wt.out
│   ├── strain_73_wt.png
│   ├── strain_73_wt_STDERR_log.out
│   ├── strain_73_wt_STDOUT_log.out
│   ├── strain_78_wt.out
│   ├── strain_78_wt.png
│   ├── strain_78_wt_STDERR_log.out
│   ├── strain_78_wt_STDOUT_log.out
│   ├── strain_93_wt.out
│   ├── strain_93_wt.png
│   ├── strain_93_wt_STDERR_log.out
│   └── strain_93_wt_STDOUT_log.out
└── metadata.xlsx
```


### 2. Preprocess Data
- Process gff/gtf file using custom script nf_core_gff_proc.sh
Bakta annotated files (tsv,gff,fna) can be found in github repository und Data/Bakta_annotated_files

Adapted from Marine Cambon: https://nfcore.slack.com/archives/CE8SSJV3N/p1679577061847429?thread_ts=1677835193.447669&cid=CE8SSJV3N

```bash
sh nf_core_gff_proc.sh -i "Data/Bakta_annotated_files/STRAINFOLDER"
```

### 3. Download Raw FASTA data from GSE277379
### 4. Use cutadapt script to trim adapter sequences 

```bash
sh cutadapt_trim.sh -i "PATHTORAWFASTQ.gz"
```

### 5. Use nextflow nf-core/RNAseq for each Strain individually:
- Create strain specific samplesheet.csv containing links to raw data downloaded from GSE277379
- Define the outdir with Strain specific subdirectories
- To run the nf-core RNA-seq pipeline (v3.14.0) using Docker, with the given sample sheet, reference genome, and specific parameters, use the following command:

```bash
nextflow run nf-core/rnaseq -r 3.14.0 -profile docker \
  --input samplesheet.csv \                                                                               # STRAIN SPECIFIC 
  --outdir nf_core_output/STRAINFOLDER \                                                                  # STRAIN SPECIFIC 
  --fasta Data/Bakta_annotated_files/STRAINFOLDER/*.fna \                                                 # STRAIN SPECIFIC 
  --gtf Data/Bakta_annotated_files/STRAINFOLDER/Processed_for_nfcoreRnaseq/filtered_transcripts.gtf \     # STRAIN SPECIFIC 
  --transcript_fasta Data/Bakta_annotated_files/STRAINFOLDER/Processed_for_nfcoreRnaseq/transcripts.fna \ # STRAIN SPECIFIC 
  --skip_umi_extract \
  --skip_trimming \
  --aligner star_salmon \
  --extra_star_align_args '--sjdbGTFfeatureExon CDS' \
  --featurecounts_feature_type transcript \
  --skip_dupradar \
  --skip_qualimap \
  --skip_rseqc \
  --featurecounts_group_type gene_id
```
### 6. Downstream Analysis and Figure Creation
Install relevant R libraries found in Libraries.R
  
### 7. Run Main.R triggering the creation of the figures found in the paper
 1. Creation of OrgDBs
 2. Creation of individual Quarto Reports, Differential Expression and Figures for each Strain
 3. Creation of Cross-Strain Quarto Report and Figures

## Data Availability:

Genomic data has been deposited in the Sequence Read Archive (SRA) and is accessible under project accession number PRJNA1158744. 
Additionally, raw RNA-seq files and TPM (transcripts per million) values of annotated genes are available through GEO under accession number GSE277379. 

## References
Harshil Patel, Phil Ewels, Alexander Peltzer, Jonathan Manning, Olga Botvinnik, Gregor Sturm, Maxime U Garcia, et al. 2024. “Nf-Core/Rnaseq: Nf-Core/Rnaseq v3.14.0 - Hassium Honey Badger.” Zenodo. https://doi.org/10.5281/ZENODO.1400710. 

James A. Fellows Yates, Jasmin Frangenberg, Anan Ibrahim, Louisa Perelo, nf-core bot, Moritz E. Beber, Hugo Tavares, Adam Talbot, Harshil Patel, and Robert Syme. 2024. “Nf-Core/Funcscan: 1.1.5 - British Beans on Toast (Patch) - 2024-03-20.” Zenodo. https://doi.org/10.5281/ZENODO.7643099. 

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim N. Artyomov, and Alexey Sergushichev. 2016. “Fast Gene Set Enrichment Analysis.” https://doi.org/10.1101/060012. 

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (12): 550. https://doi.org/10.1186/s13059-014-0550-8. 

Marc Carlson, Herve Pages. 2017. “AnnotationForge.” Bioconductor. https://doi.org/10.18129/B9.BIOC.ANNOTATIONFORGE. 

Schwengers, Oliver, Lukas Jelonek, Marius Alfred Dieckmann, Sebastian Beyvers, Jochen Blom, and Alexander Goesmann. 2021. “Bakta: Rapid and Standardized Annotation of Bacterial Genomes via Alignment-Free Sequence Identification: Find out More about Bakta, the Motivation, Challenges and Applications, Here.” Microbial Genomics 7 (11). https://doi.org/10.1099/mgen.0.000685. 

Törönen, Petri, Alan Medlar, and Liisa Holm. 2018. “PANNZER2: A Rapid Functional Annotation Web Server.” Nucleic Acids Research 46 (W1): W84–88. https://doi.org/10.1093/nar/gky350. 

Wu, Tianzhi, Erqiang Hu, Shuangbin Xu, Meijun Chen, Pingfan Guo, Zehan Dai, Tingze Feng, et al. 2021. “clusterProfiler 4.0: A Universal Enrichment Tool for Interpreting Omics Data.” The Innovation 2 (3): 100141. https://doi.org/10.1016/j.xinn.2021.100141. 



