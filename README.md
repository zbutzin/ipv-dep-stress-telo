# ipv-dep-stress-telo

### WASH Benefits Bangladesh parental IPV, depression and stress and child telomere length associations 

This is the repository for the WASH Benefits Bangladesh [NCT01590095](#https://clinicaltrials.gov/ct2/show/NCT01590095)  substudy. The primary analysis assessed the associations between maternal IPV, maternal depression, and parental stress as exposures and child telomere length as outcomes. This repository contains scripts for both the primary and post hoc analyses as well as for generating tables and figures.

### Associated protocols and datasets
The pre-specified analysis plan and the data required for the analysis will be available through the Open Science Framework: [https://osf.io/f2cm5/](https://osf.io/f2cm5/).

### WASH Benefits Package
This analysis requires [washbgam](https://github.com/washb-eed-substudies/washbgam), a package developed for WASH Benefits observational analyses. 

For all scripts, you will need to change directory statements to reflect locations of files within your local directory. Changes will also need to be made to directory statements when saving output files such as results, tables, and figures.

### Directory Structure

**`0-config`**: script for configuring data directories and loading libraries

**`src`** : analysis scripts

**`models`** : fitted models 

**`figure-data`** : data to make figures

**`figure-scripts`** : scripts to make figures from figure data

**`figures`** : resulting figures

**`results`** : analysis results (summary of fitted models)

**`table-scripts`** : scripts to make tables and resulting tables

**`tables`** : resulting tables




