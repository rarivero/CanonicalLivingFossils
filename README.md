# Canonical Living Fossils

<!-- badges: start -->
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24406663.svg)](https://doi.org/10.5281/zenodo.24406663)
[![Paper: Paleobiology (in press)](https://img.shields.io/badge/Paper-Paleobiology%20(in%20press)-red.svg)](#citation)
[![Preprint: bioRxiv](https://img.shields.io/badge/bioRxiv-preprint-black)](https://doi.org/10.1101/2024.11.14.623638)
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]
<!-- badges: end -->

## Description

This repository contains the data and scripts necessary to run the analyses contained in the paper entitled *Variable patterns of phenotypic evolution among canonical ‘living fossil’ lineages* by Rafael A. Rivero-Vega et al. (in press). The manuscript assesses whether three emblematic lineages of living fossils (coelacanths, lungfishes, and holosteans) show similar patterns of evolutionary rates over their history. Further details on the methodology can be found in the appropriate section of the paper.

The list below encompasses all of the folders, scripts, and files found within the repository and includes short explanations of each. The list broadly follows the order in which all of the data should be either created or used. Only the [R programming language](https://www.r-project.org/) is needed to replicate the analyses herein. The R scripts included automatically check whether the necessary packages are installed on the local computer and import the native R Serialized Data Files when needed (RDS files were used instead of more open formats because of a need to preserve the complex structure of the objects used and created during the analysis). Make sure the working directory is set to the top level folder containing everything in this repository.

## Contents

### Documents

- `Electronic Supplementary Material.pdf`: full electronic supplementary material for the publication in PDF format

### Folders

- `./images_*/`: contains all images used for landmarking per clade
- `./shapes_*/`: contains all ladmarking output per clade to be read into R

### Taxonomic Data

- `taxa_bodyshape.xlsx`: Microsoft Excel table with citation, clade, age, and other metadata for all taxa used in the body shape analyses
- `taxa_*_affinities.csv`: CSV with taxa and their affinity (i.e., clade/group) used in the character saturation analyses per clade
- `taxa_*_tipdates.txt`: text files with age information for all taxa per clade
- `tree_taxa_*_noshapedata.txt`: text files containing the taxa with no shape data to prune from the trees used in the shape anayses per clade
- `tree_taxa_*_outgroups.txt`: text files containing the outgroups to prune from the trees per clade

### Morphological Matrices

- `data_*_*_*.nex`: NEXUS-formatted morphological matrices for each clade 

### BEAST Input and Output

- `beast_*_morphology_*.xml`: XML files used to run BEAST per clade
- `beast_*_morphology_*_log.log`: LOG output files from BEAST per clade
- `beast_*_morphology_*.trees`: posterior trees output files from BEAST per clade
- `beast_*_morphology_*_tree-mcc.tre`: MCC consensus trees generated from trees files per clade

### R Scripts

- `0_process_raw_data.R`: pre-processing raw data step
- `1_character_saturation.R`: first character saturation step
- `2_model_fitting.R`: second model fitting step
- `3_model_bootstrapping.R`: third model bootstrapping step
- `readannotatednexus.R`: function used to read the MCC trees and preserve annotations
- `future_mvSIM.R`: custom mvSIM package that runs with futures

### Standardized Tree Data

- `tree_mcc_*_pruned.RDS`: RDS files containing the MCC trees pruned of taxa in `tree_taxa_*_outgroups.txt` for the discrete model fitting analyses per clade
- `tree_mcc_*_shapeanalyses.RDS`: RDS files containing the MCC trees pruned of taxa in `tree_taxa_*_noshapedata.txt` for the continuous model fitting analyses per clade
- `trees_posterior_*_pruned.RDS`: RDS files containing the posterior trees pruned of taxa in `tree_taxa_*_outgroups.txt` for the discrete model fitting analyses per clade
- `trees_posterior_*_shapeanalyses.RDS`: RDS files containing the prior trees pruned of taxa in `tree_taxa_*_noshapedata.txt` for the continuous model fitting analyses per clade

### Shape Data

- `data_shape_*.RDS`: RDS files containing pruned and standardized shape data per clade 

### Analytical Output

#### Character Saturation

- `saturation_*.RDS`: RDS files with calculated character saturation per clade 
- `saturation_*_vmax.RDS`: RDS files with calculated point of saturation per clade

#### Generalized Procrustes Analyses

- `gpa_*_coordinates.RDS`: RDS files with Procrustes-transformed landmark coordinate data per clade

#### Principal Component Analyses

- `pca_*_shapevalues_matrix_100.RDS`: RDS files containing PCA data representing 100% of the variance per clade

#### Model Fitting

- `fittree_*_pcamodeltests_100.RDS`: RDS files containing model fitting output using the MCC trees and PCA data representing 100% of the variance per clade
- `fittree_*_pcamodeltests_gicw_100.RDS`: RDS files containing the GIC weights for the model fitting analyses using the MCC trees and PCA data representing 100% of the variance per clade
- `fittrees_*_pcamodeltests_100.RDS`: RDS files containing model fitting output using the posterior trees and PCA data representing 100% of the variance per clade
- `fittrees_*_pcamodeltests_gicw_100.RDS`: RDS files containing the GIC weights for the model fitting analyses using the posterior trees and PCA data representing 100% of the variance per clade

#### Bootstrapping

- `fitboot_*_pcamodeltests_100.RDS`: RDS files containing the bootstrapping results using the model fitting output, MCC trees, and PCA data representing 100% of the variance per clade
- `fitboot_*_pcamodeltests_stats_100.RDS`: RDS files containing statistics calculated from the bootstrapping results using the model fitting output, MCC trees, and PCA data representing 100% of the variance per clade

## Citations

Rafael A. Rivero-Vega, Jacob S. Berv, John T. Clarke, Matt Friedman. *Inconsistent patterns of phenotypic evolution among canonical ‘living fossil’ lineages*. **Paleobiology**. In press. [Preprint on bioRxiv](https://doi.org/10.1101/2024.11.14.623638)

Rafael A. Rivero-Vega, 2024. *From Living Fossils to Adaptive Radiations: Examining Phenotypic Evolution in Lobe-Finned Fishes*. [Doctoral Dissertation](https://dx.doi.org/10.7302/24028)

Rafael A. Rivero-Vega. 2024. *From Living Fossils to Adaptive Radiations: Examining Phenotypic Evolution in Lobe-Finned Fishes* \[Data set\]. **University of Michigan - Deep Blue Data**. [Dissertation Dataset](https://doi.org/10.7302/4mcp-xc40)

## License

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
