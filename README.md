# Genomic and Transcriptomic Adaptation to Chlorhexidine in Streptococcus spp.

Bernd Daller1, David L. Auer2, Wolfgang Buchalla2, Sibylle Bartsch3, André Gessner1, Nicholas S. Jakubovics4, Ali Al-Ahmad3, Andreas Hiergeist1,#, Fabian Cieplik3,#,\*

1 Institute of Clinical Microbiology and Hygiene, University Hospital Regensburg, Regensburg, Germany

2 Department of Conservative Dentistry and Periodontology, University Hospital Regensburg, Regensburg, Germany

3 Department of Operative Dentistry and Periodontology, Center for Dental Medicine, Medical Center – University of Freiburg, Medical Faculty, University of Freiburg, Freiburg i. Br., Germany

4 School of Dental Sciences, Faculty of Medical Sciences, Newcastle University, Newcastle upon Tyne, UK

------------------------------------------------------------------------

This study investigated CHX adaptation in three Streptococcus spp. clinical isolates and their CHX-adapted counterparts using whole genome sequencing and transcriptomic profiling.

## Follow these steps to replicate the RNAseq Analysis:

### 1. Clone repository containing relevant Data and scripts:

``` bash
├── 01_Main_download_and_process.sh
├── 02_Main_data_analysis.R
├── bash
│   ├── download_rnaseq_data.sh
│   ├── preprocess_gff.sh
│   ├── run_nfcore_rnaseq.sh
│   ├── trim_adapters.sh
├── Data
│   ├── Antibiotic_Resistance
│   ├── Bakta_annotated_files
│   ├── GO_Check
│   ├── GO_Pannzer2
├── OrgDbs
│   ├── org.SmitisS93wt.eg.db
│   ├── org.SsalivariusS73wt.eg.db
│   └── org.SvestibularisS78wt.eg.db
├── R
│   ├── CreateOrgDb.R
│   ├── CrossStrain_Analyis.qmd
│   ├── CrossStrain_helper_functions.R
│   ├── DifferentialExpression_helper_functions.R
│   ├── DifferentialExpression.qmd
│   ├── Libraries.R
│   └── _quarto.yml
├── README.md
├── StrepCHX2024.Rproj
```

### 2. Preprocess Data using `01_Main_download_and_process.sh`

-   `cd` to the location where you cloned the git repository
-   Execute `bash 01_Main_download_and_process.sh`

### 3. Downstream Analysis using `02_Main_data_analysis.R`

-   Creation of OrgDBs
-   Creation of individual Quarto Reports, Differential Expression Analysis and Figures for each Strain
-   Creation of Cross-Strain Quarto Report and Figures

#### Structure after execution of 01_Main_download_and_process.sh and 02_Main_data_analysis.R:

``` bash
├── 01_Main_download_and_process.sh
├── 02_Main_data_analysis.R
├── bash
│   ├── download_rnaseq_data.sh
│   ├── preprocess_gff.sh
│   ├── run_nfcore_rnaseq.sh
│   ├── trim_adapters.sh
├── Data
│   ├── Antibiotic_Resistance
│   ├── Bakta_annotated_files
│   ├── GO_Check
│   ├── GO_Pannzer2
│   └── PRJNA1162077_data*
├── nf_core_output*
│   ├── Strain_73*
│   ├── Strain_78*
│   └── Strain_93*
├── OrgDbs
│   ├── org.SmitisS93wt.eg.db
│   ├── org.SsalivariusS73wt.eg.db
│   └── org.SvestibularisS78wt.eg.db
├── output*
│   ├── CrossStrain_Analyis*
│   ├── Quarto_Reports*
│   ├── Strain_73*
│   ├── Strain_78*
│   └── Strain_93*
├── R
│   ├── CreateOrgDb.R
│   ├── CrossStrain_Analyis.qmd
│   ├── CrossStrain_helper_functions.R
│   ├── DifferentialExpression_helper_functions.R
│   ├── DifferentialExpression.qmd
│   ├── Libraries.R
│   └── _quarto.yml
├── README.md
├── StrepCHX2024.Rproj
├── work_strain_73*
├── work_strain_78*
├── work_strain_93*
```

## Data Availability:

Genomic data has been deposited in the Sequence Read Archive (SRA) and is accessible under project accession number PRJNA1158744. Additionally, raw RNA-seq files and TPM (transcripts per million) values of annotated genes are available through GEO under accession number GSE277379.\
Nevertheless the script should download and process the data automatically.

## References

Harshil Patel, Phil Ewels, Alexander Peltzer, Jonathan Manning, Olga Botvinnik, Gregor Sturm, Maxime U Garcia, et al. 2024. “Nf-Core/Rnaseq: Nf-Core/Rnaseq v3.14.0 - Hassium Honey Badger.” Zenodo. <https://doi.org/10.5281/ZENODO.1400710>.

James A. Fellows Yates, Jasmin Frangenberg, Anan Ibrahim, Louisa Perelo, nf-core bot, Moritz E. Beber, Hugo Tavares, Adam Talbot, Harshil Patel, and Robert Syme. 2024. “Nf-Core/Funcscan: 1.1.5 - British Beans on Toast (Patch) - 2024-03-20.” Zenodo. <https://doi.org/10.5281/ZENODO.7643099>.

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim N. Artyomov, and Alexey Sergushichev. 2016. “Fast Gene Set Enrichment Analysis.” <https://doi.org/10.1101/060012>.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (12): 550. <https://doi.org/10.1186/s13059-014-0550-8>.

Marc Carlson, Herve Pages. 2017. “AnnotationForge.” Bioconductor. <https://doi.org/10.18129/B9.BIOC.ANNOTATIONFORGE>.

Schwengers, Oliver, Lukas Jelonek, Marius Alfred Dieckmann, Sebastian Beyvers, Jochen Blom, and Alexander Goesmann. 2021. “Bakta: Rapid and Standardized Annotation of Bacterial Genomes via Alignment-Free Sequence Identification: Find out More about Bakta, the Motivation, Challenges and Applications, Here.” Microbial Genomics 7 (11). <https://doi.org/10.1099/mgen.0.000685>.

Törönen, Petri, Alan Medlar, and Liisa Holm. 2018. “PANNZER2: A Rapid Functional Annotation Web Server.” Nucleic Acids Research 46 (W1): W84–88. <https://doi.org/10.1093/nar/gky350>.

Wu, Tianzhi, Erqiang Hu, Shuangbin Xu, Meijun Chen, Pingfan Guo, Zehan Dai, Tingze Feng, et al. 2021. “clusterProfiler 4.0: A Universal Enrichment Tool for Interpreting Omics Data.” The Innovation 2 (3): 100141. <https://doi.org/10.1016/j.xinn.2021.100141>.
