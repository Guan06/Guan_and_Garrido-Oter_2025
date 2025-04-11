## Novel microbial diversity and network analyses identify drivers of dynamics in plant microbiota 

### Abstract
Microbial communities are key components of most ecosystems, where interactions among microbes shape fundamental processes that support biodiversity and productivity. With the rapid development of sequencing technologies, an increasing number of microbiome datasets have been generated. However, standard analyses of these datasets often compare community composition between samples but overlook how microbes co-occur and potentially influence each other. To address this problem, we developed a generalizable computational framework that integrates compositional and co-occurrence network analyses, which we then applied to extensive microbial amplicon datasets. Here, we focus on the plant microbiota, which typically exhibits high diversity and remains challenging to characterize due to the large number of low-abundance taxa. We show that identifying a subset of representative microbial taxa captures the overall community structure while increasing the statistical power of downstream analyses. From these representative taxa, we inferred a large-scale co-occurrence network, from which microbes with co-varying abundances were clustered into units for diversity measurement. This approach not only reduces unexplained variance in diversity assessments but also captures the key microbe–microbe relationships that govern assembly patterns. Furthermore, we introduced a bootstrap- and permutation-based statistical approach to compare microbial networks from diverse conditions or environments. We demonstrate that our method robustly distinguishes meaningful differences and pinpoints the specific microbes and microbial features that drive those differences. These results highlight the importance of incorporating microbe-microbe interactions in microbiota studies, leading to more accurate and ecologically meaningful insights than traditional compositional comparisons alone. Our framework is made freely available as an R package, ‘mina’, enabling researchers to analyze microbial network architecture, identify condition-specific microbial interactions, and ultimately gain a deeper understanding of community ecology. With its broad applicability beyond plant systems, this package provides a valuable tool for leveraging microbiome data across disciplines, from agricultural settings to ecosystem resilience and human health.

### Data and code availability
All the data used in this study is published and cited. Scripts used to generate the analyses and figures could be found in this repository. The R package ‘mina’ can be downloaded and installed from either Bioconductor [10.18129/B9.bioc.mina](https://www.bioconductor.org/packages/release/bioc/html/mina.html) or GitHub repository [mina](https://github.com/Guan06/mina).

### Repository structure

In this repository, you will find the following subfolders:

- __00.data__
	This folder contains metadata of samples, relative abundance of bacteria and fungi in ASV table format processed with the DADA2-based pipeline [here](https://github.com/Guan06/DADA2_pipeline), taxonomic assignments of ASVs, and calculated alpha diversity (average of 999 bootstrap).

- __01.common_scripts__
	This folder holds functions and plot settings used in each subfolder to reproduce analyses and figures.

- __02.figure_1_and_S1_to_S4__
	This folder contains scripts to process and visualize results for 
        1) Main Figure 1.
		2) Supplementary Figures 1-4. 
		3) Statistical results and script output are also included.

Additional subfolders named in the same way contains similar files with a consistent structure.

### Authors and contributors 

This repository is written by [Rui Guan](https://github.com/Guan06) with contributions from [Ruben Garrido-Oter](https://github.com/garridoo).

Follow me on Bluesky at [guan06rui.bsky.social](https://bsky.app/profile/guan06rui.bsky.social).
