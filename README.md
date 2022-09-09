# Targeted Metabolomics Processing Pipeline

The Targeted Pipeline is a process for importing, cleaning, QCing, BMISing, and quantifying (where possible) targeted data produced by the Ingalls Lab.

It is run by a "Control" script that moves through the pipeline step by step, accounting for mass spectrometer instrument type (TQS/QE), run type (HILIC/Reverse phase), ion mode (positive/negative), and metabolomics processing software (Skyline or MsDial).

The Control.Rmd script contains detailed instructions on each section, and *the user must read carefully through each section while progressing through the pipeline*.

The pipeline contains four major sections:

1.  Import and cleaning/rearranging of data.
2.  Quality control using user-defined parameters.
3.  Applying Best-Matched Internal Standard (B-MIS).
4.  Quantifying peak area to umol/vial when possible.

This code should be viewed as a template only. It will need to be adjusted for each run, and this repository exists only as a training tool!

## Usage

The following packages are required for the pipeline.

``` r
library(anytime)
library(rlist)
library(tidyr)
library(tidyverse)
```

The first few chunks of the Control.Rmd script will create empty folders for different types of data (raw, processed, extra, intermediate). Some of this data will be produced by the script (processed, intermediate) and some needs to be moved into those folders before proceeding. Because this is a template script, the data has been modified for easy input, which will NEVER be the case in a real run.

All the necessary example data can be downloaded here: <https://drive.google.com/drive/folders/1E-B4dyDTkwOfupHJK0uL8BgKj586OU31>

## Visualization

![Image Title](visual/Targeted_Pipeline_Visualization.pdf)

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
