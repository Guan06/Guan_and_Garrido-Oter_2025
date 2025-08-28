## Integrated diversity and network analyses reveal drivers of microbiome dynamics

### Abstract
Microbial communities are key components of ecosystems, where interactions among microbes drive biodiversity and productivity. Increased number of microbiome datasets are available owing to advance in sequencing, yet standard analyses often focus on community composition, neglecting the complex interactions between co-occurring microbes. To address this, we developed a computational framework integrating compositional and co-occurrence network analyses. We applied to extensive microbial amplicon datasets, focusing on plant microbiota, which typically exhibits high diversity and remains challenging to characterize due to the large number of low-abundance taxa. We show that identifying a subset of representative microbial taxa captures the overall community structure and increases the statistical power. From these taxa, we inferred a large-scale co-occurrence network and clustered microbes with co-varying abundances into units for diversity measurement. This approach not only reduces unexplained variance in diversity assessments but also captures the key microbe–microbe relationships that govern assembly patterns. Furthermore, we introduced a bootstrap- and permutation-based statistical approach to compare microbial networks from diverse conditions. Our method robustly distinguishes meaningful differences and pinpoints specific microbes and features driving those differences. These results highlight the importance of incorporating microbe-microbe interactions in microbiota studies, leading to more accurate and ecologically meaningful insights. Our framework, available as an R package ('mina'), enables researchers to identify condition-specific interactions via network comparison, and gain a deeper understanding of community ecology. With broad applicability beyond plant systems, this package provides a valuable tool for leveraging microbiome data across disciplines, from agriculture to ecosystem resilience and human health.

### Importance
Understanding microbiome dynamics requires capturing not only changes in microbial composition but also interactions between community members. Traditional approaches frequently overlook microbe-microbe interactions, limiting their ecological interpretation. Here, we introduce a novel computational framework that integrates compositional data with network-based analyses, significantly improving detection of biologically meaningful patterns in community variation. Application of this framework to a large dataset from the plant microbiota, we identify representative groups of interacting microbes driving differences across microhabitats and environmental conditions. Our analysis framework, implemented in an R package ‘mina’, provides robust tools allowing researchers to assess statistical differences between microbial networks and detect condition-specific interactions. Broadly applicable to microbiome datasets, we aim at enabling advances in our understanding of microbial interactions within complex communities.

### Data and code availability
All the data used in this study is published and cited. Scripts used to generate the analyses and figures could be found in this repository. The R package ‘mina’ can be downloaded and installed from either Bioconductor [10.18129/B9.bioc.mina](https://www.bioconductor.org/packages/release/bioc/html/mina.html) or GitHub repository [mina](https://github.com/Guan06/mina).

### Repository structure

In this repository, you will find the following subfolders:

- __00.data__
	This folder contains metadata of samples, relative abundance of bacteria and fungi in ASV table format processed with the DADA2-based pipeline [here](https://github.com/Guan06/DADA2_pipeline), taxonomic assignments of ASVs, and calculated alpha diversity (average of 999 bootstrap).

- __01.common_scripts__
	This folder holds functions and plot settings used in each subfolder to reproduce analyses and figures.

- __02.figure_1_and_S1_to_S7__
	This folder contains scripts to process and visualize results for 
        1) Main Figure 1.
		2) Supplementary Figures 1-7. 
		3) Statistical results and script output are also included.

Additional subfolders named in the same way contains similar files with a consistent structure.

### Authors and contributors 

This repository is written by [Rui Guan](https://github.com/Guan06) with contributions from [Ruben Garrido-Oter](https://github.com/garridoo).

Follow me on Bluesky at [guan06rui.bsky.social](https://bsky.app/profile/guan06rui.bsky.social).

### Citation

to be added..

### Acknowledgements

Acknowledgements
We would like to thank Prof. Dr. Paul Schulze-Lefert and Prof. Dr. Gunnar W. Klau for their advice throughout the duration of this project. We thank [Dr. Yulong Niu](https://github.com/YulongNiu) for helping with submitting the R package ‘mina’.
