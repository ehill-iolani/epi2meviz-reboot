# EPI2MEViz - Reboot

A R based shiny app that is intended to be used as a companion analysis tool for both wf-16s and wf-metagenomics EPI2ME workflows from Oxford Nanopore Technologies.

This is the sucessor to the orignal EPI2MEViz released in 2023 prior to the large EPI2ME update from Oxford Nanopore Technologies.

This app will take the .csv file downloaded from the EPI2ME analysis and will conduct 4 analyses:
1) Generate rarefaction curves to determine if the run has achieved the proper sampling depth
2) Calculate the relative abundance of each taxa detected in the EPI2ME results and generate the corresponding plot. The plot will allow you to view relative abundance at different taxonomic levels
3) Calculate alpha diversity using the Simpson's index to compare the diversity between samples.
4) Conduct a Bray Curtis dissimilarity principle coordinate analysis (PCoA) and plot PCoA 1 vs PCoA 2 (only works if there are > 2 barcodes in the analysis)

You will need to upload the .csv containing classifications file downloaded from either wf-16s or wf-metagenomics analysis, and optionally, a medata files that corresponds with the barcodes in the analysis.
