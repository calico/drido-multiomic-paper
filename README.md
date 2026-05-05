# drido-multiomic-paper

_Companion R code for "A multiomic lifespan signature in genetically diverse, diet-restricted mice"_

### Instructions 

This package contains the code necessary to perform the statistical analysis and
generate the figures shown in the associated publication. User may need to access
raw data files from online repositories to run complete analyses. This code does
not build as a package but instead sources function and data files as listed below.
Recommended use is via forking and cloning:

With SSH: `git clone git@github.com:calico/drido-multiomic-paper.git`
With PAT: `git clone https://<GITHUB_USERNAME>:<PAT>@github.com/calico/drido-multiomic-paper.git`

## Contents
  
- **/data:**
  - Contains files required for lipid normalization (.msp files)
  name keys for GO enrichment, etc.
  
- **/R:**
  - R script files to functions sourced in `/inst` files
  
- **/inst:**
  - **/extdata:** R data objects (`.Rds`) generated from statistical analysis
  scripts and used to make figures
  - **/figures:** R markdown (`.Rmd`) files used to generate figures in publication
  - **/scripts:** R scripts used for normalization of metabolomics and lipidomics
  data, statistical analyses, and parsing of data. 
    - Proteomics data was processed using [msTrawler](https://github.com/calico/msTrawler)
    - Lipidomics data processing herein requires [mzkitcpp](https://github.com/calico/mzkitcpp)
    - Metabolomics data was annotated with [MAVEN](https://github.com/eugenemel/maven)
    - Other normalization processing herein requires [claman](https://github.com/calico/claman) and [romic](https://github.com/calico/romic)
  - **/supp_tables:** Subset of published supplemental data used for making figures
  or deploying shiny app
  
- **/shiny-app:**
  - Files for local deployment of interactive data explorer app
  - To deploy from the repository root, run `shiny::runApp('shiny-app')`
  

