Targeted Metabolomics Processing Pipeline
================

> Microbial Metabolomics Research Center

[![PRs
welcome](https://img.shields.io/badge/PRs-welcome-ff69b4.svg)](https://github.com/nhn/tui.editor/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22)

## Table of Contents

-   [Description](#-Description)
-   [Packages and External Processes](#-Packages-and-External-Processes)
-   [Usage](#-Usage)
-   [Visualization](#-Visualization)

------------------------------------------------------------------------

### Description

**Important:** *This repository is a template and training tool! Each
dataset is different; while this pipeline handles many of the scenarios
that arise within the data wrangling process, the user will need to
adjust the parameters for their specific input data.*

The Targeted Pipeline is a combined process for importing, tidying,
Quality Controlling, BMISing (see below for more details on this step),
and quantifying (where possible) targeted metabolomics data produced by
the Ingalls Lab.

It is run by a **“Control.Rmd”** markdown script that moves through the
pipeline step by step, accounting for two mass spectrometer instrument
types
([TQS](https://www.waters.com/waters/en_US/Xevo-TQ-S/nav.htm?cid=10160596&locale=en_US)
or
[QE](https://www.thermofisher.com/order/catalog/product/IQLAAEGAAPFALGMBDK)),
run type (HILIC/Reverse phase), ion mode (positive/negative), and
metabolomics processing software
([MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html) or
[Skyline](https://skyline.ms/wiki/home/software/Skyline/page.view?name=tutorial_hi_res_metabolomics)).

The Control.Rmd script contains detailed instructions, and *the user
must read carefully through each section description while progressing
through the pipeline*.

The pipeline contains four major sections:

1.  Import and cleaning/rearranging of data.
2.  Quality control using user-defined parameters.
3.  Applying a Best-Matched Internal Standard (B-MIS).
4.  Quantifying peak area to umol/vial when possible.

### 📦 Packages and External Processes

Below is a list of packages and external sofwares that this pipeline
utilizes. Please ensure you have the package(s) installed and have
access to the tools listed below.

| Name                                                                    | Description                                                                                                                                                                                                                                                                                                                                                                                |
|:------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`tidyverse`](https://www.tidyverse.org/)                               | Data science package in R.                                                                                                                                                                                                                                                                                                                                                                 |
| [`anytime`](https://cran.r-project.org/web/packages/anytime/index.html) | General purpose date converter. This package helps the user track analysis reruns and produce unique identifiers for files.                                                                                                                                                                                                                                                                |
| [`rlist`](https://cran.r-project.org/web/packages/rlist/index.html)     | Non-tabular data manipulation toolbox.                                                                                                                                                                                                                                                                                                                                                     |
| [`BMIS`](https://github.com/IngallsLabUW/B-MIS-normalization)           | R script for batch-specific normalization of raw peak areas. Matches measured metabolites with isotope-labeled internal standards that behave similarly during the analysis, nicknamed Best-Matched Internal Standard or B-MIS. This pipeline contains a modified script taken from the original repository. Published paper here: <https://pubs.acs.org/doi/10.1021/acs.analchem.7b04400> |

### Usage

The first few chunks of the Control.Rmd script will create empty folders
for different types of data (raw, processed, extra, intermediate). Some
of this data will be produced by the script (processed, intermediate)
and some needs to be moved into those folders before proceeding. Because
this is a template script, the data has been modified for easy input,
which will NEVER be the case in a real run.

All the necessary example data can be downloaded here:
<https://drive.google.com/drive/folders/1E-B4dyDTkwOfupHJK0uL8BgKj586OU31>

## Visualization

![Click on me to see a visual layout of the Targeted
Pipeline!](visual/Targeted_Pipeline_Visualization.pdf)

------------------------------------------------------------------------

## 🔧 Pull Requests

Pull requests are welcome. For major changes, please open an issue
first.

## 💬 Contributors

-   [Regina Lionheart](https://github.com/R-Lionheart)
-   [Katherine Heal](https://github.com/kheal)
-   [Angie Boysen](https://github.com/Angie-B)
