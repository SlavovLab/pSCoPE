---
layout: default
title: Download data
nav_order: 2
permalink: mass-spec/data
description: "download single-cell proteomics data from pSCoPE, a prioritized SCoPE-MS method"
---
{% include social-media-links.html %}

<!--# Download single-cell protein and RNA data-->

&nbsp;

[Processed pSCoPE Data]({{site.baseurl}}#processed-single-cell-protein-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[RAW pSCoPE Data]({{site.baseurl}}#RAW-single-cell-protein-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
<!--[10x Genomics Data]({{site.baseurl}}#single-cell-RNA-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }-->

&nbsp;

<h2 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >pSCoPE data processed to ASCII text matrices</h2>
<h3 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >Benchmarking experiments: Figure 1b/e data</h3>

* [Peptides-raw.csv](https://drive.google.com/file/d/1gbbn_SgwaPYAYG3iiQqpEHSw34hpKp5M/view?usp=sharing)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/). 
  <!-- See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details. -->  

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1l1DvtdYrfj6rjq1TZKcsWqG4iYBOtpKq/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/1J7sBHYwo5669w8igY4ldojPcLwuqX7I-/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each row corresponds to a single cell annotated with relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;

<!--* [sdrf_meta_data.tsv](https://drive.google.com/file/d/1T8BTfNDlYQkBTs8La6YRSCyD1RwNTvqk/view?usp=sharing)
   -  Meta data following the [Sample to Data file format (SDRF) for Proteomics project guidelines](https://github.com/bigbio/proteomics-metadata-standard) for  for all single cells used in analysis constituting all figures.-->

&nbsp;

<!--* [Joint protein-RNA data](https://drive.google.com/file/d/130FWc-s-Pd-mx3ymg22bI1qH5fiT7Ktv/view?usp=sharing)
   - `Gene` **x**  `single cells`. Both sets imputed and batch-corrected separately then combined, taking only genes common to both data sets. Uniprot accession numbers used to denote gene.-->

<!--&nbsp;-->

<!--* [Signal-to-noise data](https://drive.google.com/file/d/16dmI7qNdpJlPOn83dOZFhHfXv0du5Dip/view?usp=sharing)
  - `Peptides` and `Proteins` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. The quantitation is the Signal-to-noise (S/N) ratio for each single cell's corresponding reporter ion extracted from the RAW file. The single cell identification numbers are [mapped](https://drive.google.com/file/d/1PUfiGhmInYP3JW5Xoul7Tikl9RSHyQcN/view?usp=sharing) to cell type and RAW file. Complete extracted S/N for each RAW file can be found [here](https://drive.google.com/drive/folders/18_BQ15_JQKzbDt1JZo36MaJuOhN3tJCX?usp=sharing).  -->

<h3 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >Benchmarking experiments: Figure 1c/d data</h3>

* [Peptides-raw.csv](https://drive.google.com/file/d/1VS-ko7rDsy0t2V5JNfB049txTLpzGPr0/view?usp=sharing)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/). 
  <!-- See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details. -->

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1PSg-FA5eqYeLEekVZ8Yr_nx902K4Aink/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/1z6vah8XfzqxEM62oFMmGP7q7UDYhfS6D/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each row corresponds to a single cell annotated with relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;

<!-- * [sdrf_meta_data.tsv](https://drive.google.com/file/d/1T8BTfNDlYQkBTs8La6YRSCyD1RwNTvqk/view?usp=sharing)
   -  Meta data following the [Sample to Data file format (SDRF) for Proteomics project guidelines](https://github.com/bigbio/proteomics-metadata-standard) for  for all single cells used in analysis constituting all figures. -->  

<h3 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >Bone-marrow-derived macrophage experiments: Figures 2, 3, and 4</h3>

* [Peptides-raw.csv](https://drive.google.com/file/d/1mCNeDcxUT5eWKwSruC-aThnKM9PNtEcO/view?usp=sharing)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. 
  <!-- See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details. -->  

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1qHU5wtXoKxBcZ73QPWSUavOeO2xk1zvR/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/1OHEf8PQ7REerh0kFlnxj3bHFQ2HTmRKL/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each row corresponds to a single cell annotated with relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;

<!-- * [sdrf_meta_data.tsv](https://drive.google.com/file/d/1T8BTfNDlYQkBTs8La6YRSCyD1RwNTvqk/view?usp=sharing)
   -  Meta data following the [Sample to Data file format (SDRF) for Proteomics project guidelines](https://github.com/bigbio/proteomics-metadata-standard) for  for all single cells used in analysis constituting all figures. -->


&nbsp;
<!-- * [DART-ID input](https://drive.google.com/drive/folders/1ohLco5KHX95jyXIZUAZDvrrbip1RzZ_1?usp=sharing) -->


&nbsp;
<!-- * [GSEA: GOrilla output](https://drive.google.com/drive/folders/1DCp_euY0Cj_NWWG5xQsx7CTN3ju5LI_O?usp=sharing)
&nbsp; -->

<!-- * [Minimal data files](https://drive.google.com/drive/folders/10pOMMlxHsFIyPa9X2auq6xKJssqFgo-D?usp=sharing) necessary for generating Peptides-raw.csv and Proteins-processed.csv -->

<h2 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >Files needed to regenerate preprint figures</h2>

&nbsp;

* [Additional data files](https://drive.google.com/file/d/1WRVpDCJPxYjfOX-QWvoq86RqfKbmVfGc/view?usp=sharing) necessary for generating figures from the pSCoPE preprint. 
<!-- [SCoPE2 article](https://doi.org/10.1101/665307). -->

<!--&nbsp;-->

<!-- * [Processed Data](https://drive.google.com/drive/folders/1NJODxiKrnfW2_nTP-_n_UDvIpwcDEz4C?usp=sharing) from the [second version (v2) of the SCoPE2 preprint](https://www.biorxiv.org/content/10.1101/665307v3) -->

<!--&nbsp;-->

<!-- * [Processed Data](https://drive.google.com/open?id=1cMQ-SIGpHwSfx9wJF2fIa-t8yX329LPM) from the [first version (v1) of the SCoPE2 preprint](https://www.biorxiv.org/content/10.1101/665307v1) -->

<!--&nbsp;-->

<!-- * [Single cell proteomics data processing](https://uclouvain-cbio.github.io/scp/): The analysis of the data described here has been replicated by Christophe Vanderaa and Laurent Gatto with the scp Bioconductor package: The scp package is used to process and analyze mass spectrometry-based single cell proteomics data and is freely available from their [Github repository](https://github.com/UCLouvain-CBIO/scp/). The scp package and the 
replication are described in this [video](https://youtu.be/XMxZkw8yorY). -->

&nbsp;

<!-- ## DO-MS reports of SCoPE2 data

To facilitate the exploration of the SCoPE2 data, we plotted the distributions of important and informative features of the LC-MS/MS data using the methodology of [Data-driven Optimization of MS (DO-MS)](https://do-ms.slavovlab.net) developed by Huffman et al., _J. of Proteome Research_, 2019, DOI: [10.1021/acs.jproteome.9b00039](https://doi.org/10.1021/acs.jproteome.9b00039) -->


<!-- * [Set 1]({{site.baseurl}}A1_glance/index.html) -->
<!-- * [Set 1]({{site.baseurl}}B1_glance/index.html) -->
<!-- * [Set 2]({{site.baseurl}}B2_glance/index.html) -->



&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="RAW-single-cell-protein-data" >pSCoPE RAW data and search results from MaxQuant</h2>
The Massive repository below contains RAW mass-spectrometry data files generated by a Q exactive instrument as well as the search results from analyzing the RAW files by [MaxQuant](https://www.maxquant.org/)  and by [DART-ID](https://dart-id.slavovlab.net/). Please consult the linked files for more information on the included files:

<!--* [MaxQuant results descriptions](https://drive.google.com/open?id=1qXThKpGPx1tBcxvYFvNM0zCSeyILDzE6) -->

* [Raw file descriptions](https://docs.google.com/spreadsheets/d/1EoBFPIgXXSYqP5khTAp7LMG75HYGK2uT/edit?usp=sharing&ouid=109814487119977139380&rtpof=true&sd=true)

* **MassIVE Repository:**
  - [**http:**  MSV000089055](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b15cafc7489147e99b93bd7c718388b2)
  - [**ftp:** &nbsp; MSV000089055](ftp://massive.ucsd.edu/MSV000089055/)

<!-- * **MassIVE Repository 2:**
  - [**http:**  MSV000084660](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000084660)
  - [**ftp:** &nbsp; MSV000084660](ftp://massive.ucsd.edu/MSV000084660) -->



&nbsp;


<!--<h2 style="letter-spacing: 2px; font-size: 26px;" id="single-cell-RNA-data" >scRNA-seq 10x Genomics RAW and processed data</h2>

A cellular mixture identical to that used for the single-cell proteomics was assessed with scRNA-seq using 10x Genomics Chromium platform and the Single Cell 3’ Library & Gel Bead Kit (v2). Two biological replicates of the cell suspension (stock concentration: 1200 cells/μl) were loaded into independent lanes of the device. An average number of about 10,000 cells/lane were recovered. Following the library preparation and sample QC by Agilent BioAnalyzer High Sensitivity chip, the two libraries were pooled together, quantified by KAPA Library Quantification kit and sequenced using the Illumina Novaseq 6000 system (Nova S1 100 flow cell) with the following run parameters: Read1: 26 cycles, i7 index: 8 cycles, Read2: 93 cycles.

* **GEO Repository**
  - [**http:**  GSE142392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142392)

  &nbsp;  

* **Processed data matrices: `Transcript` x  `single cells` in UMI counts.**
  - [mRNA biological replicate one](https://drive.google.com/open?id=1cN6UgSrZfqKdOjwJ0VyEYp6m_Fy9eANR)

  - [mRNA biological replicate two](https://drive.google.com/open?id=1cuoYiqKgzVnUoFnFmrpXWKVfaFwiboeo) -->

&nbsp;  

&nbsp;

&nbsp;  

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

