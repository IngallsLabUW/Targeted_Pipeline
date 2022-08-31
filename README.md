# Targeted Metabolomics Processing Pipeline

### aka, the Targeted Pipeline

The Targeted Pipeline is a process for importing, cleaning, QCing, BMISing, and quantifying (where possible) targeted data produced by the Ingalls Lab.

It is run by a "Control" script that moves through the pipeline step by step, accounting for mass spectrometer instrument type (TQS/QE), run type (HILIC/Reverse phase), ion mode (positive/negative), and metabolomics processing software (Skyline or MsDial).

The Control.Rmd script contains detailed instructions on each section, and *the user must read carefully through each section while progressing through the pipeline*.

The four major sections are: Section I: Import and cleaning/rearranging of data. Section II: Quality control using user-defined parameters. Section III: Applying Best-Matched Internal Standard (B-MIS). Section IV: Quantifying peak area to umol/vial when possible.

## Usage
The following packages are required for the pipeline. 

``` bash
library(anytime)
library(rlist)
library(tidyr)
library(tidyverse)
```

## Usage

``` python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
