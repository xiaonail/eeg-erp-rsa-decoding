# EEG N170 ERP/RSA Pipeline

MATLAB-based EEG analysis pipeline for N170 ERP analysis and representational similarity analysis (RSA).

The pipeline includes dataset discovery, preprocessing, ERP analysis, time-frequency analysis, RSA, and summary reporting. 

## Dataset

This pipeline is based on the **ERP CORE N170 face-perception dataset** from the ERP CORE project.

The present pipeline uses the **N170 face-perception task**, in which participants viewed:

- faces;
- cars;
- scrambled faces;
- scrambled cars.

The main ERP contrast in this repository focuses on **faces vs. cars** in the N170 time window. The RSA analysis compares EEG dissimilarity structure with a simple face/car category model.

Raw EEG data are not redistributed in this repository. Please download the dataset from the official ERP CORE/OSF source and follow the original data-use terms.

## Repository Structure

```text
code/        Main analysis scripts
functions/   Helper functions used by the pipeline
results/     Example output figures
LICENSE      Repository license
```

## Requirements

- MATLAB
- EEGLAB
- EEGLAB ICLabel plugin
- EEGLAB clean_rawdata plugin
- FieldTrip, required only for the optional time-frequency analysis

The required MATLAB toolboxes should be installed separately and added to the MATLAB path before running the pipeline.

## Setup

Copy the configuration template:

```text
code/config_template.m
```

Rename it to:

```text
code/config.m
```

Then edit `config.m` to define the local paths for the raw data, processed data, results folder, EEGLAB, and FieldTrip.

## Running the Pipeline

Run the scripts in order from the `code/` folder:

```matlab
step0_setup
step1_discover
step2_preprocess
step3_erp
step4_timefreq
step5_rsa
step6_summary
```

## Outputs

The pipeline generates ERP plots, scalp topographies, time-frequency results, RSA figures, and a summary of key analysis outputs.

## Data and Reproducibility

Raw EEG data, processed EEG files, and generated MATLAB result files should not be committed to GitHub. The repository is intended to store analysis code, documentation, and lightweight example figures only.

## Citation

If using the ERP CORE dataset, cite the original ERP CORE paper:

Kappenman, E. S., Farrens, J. L., Zhang, W., Stewart, A. X., & Luck, S. J. (2021). *ERP CORE: An open resource for human event-related potential research*. NeuroImage, 225, 117465. https://doi.org/10.1016/j.neuroimage.2020.117465

## Notes

This code is intended for reproducible EEG research analysis. Results should be interpreted with attention to sample size, preprocessing choices, statistical correction, and dataset-specific assumptions.

## License

The code is released under the MIT License.
