# EEG N170 ERP RSA decoding Pipeline

MATLAB and Python EEG analysis pipeline for N170 ERP analysis, representational similarity analysis, multivariate decoding, and spectral parameterisation.

The pipeline covers dataset discovery, preprocessing, ERP analysis, time-frequency analysis, RSA, decoding, periodic and aperiodic spectral components, and summary reporting.

## Dataset

This pipeline is based on the **ERP CORE N170 face-perception dataset** from the ERP CORE project.

The present pipeline uses the **N170 face-perception task**, in which participants viewed:

- faces;
- cars;
- scrambled faces;
- scrambled cars.

The main ERP contrast in this repository focuses on **faces vs. cars** in the N170 time window. The RSA analysis compares EEG dissimilarity structure with a simple face/car category model. The decoding analysis classifies single face and car trials at each time point. The spectral parameterisation separates the aperiodic exponent and offset of the post-stimulus power spectrum from oscillatory peaks and compares the exponent between faces and cars.

## Requirements

MATLAB steps:

- MATLAB
- EEGLAB
- FieldTrip

The required MATLAB toolboxes should be installed separately and added to the MATLAB path before running the pipeline.

Python decoding and spectral parameterisation steps:

- Python 3.10 or later
- MNE-Python, scikit-learn, SciPy, NumPy, Matplotlib, FOOOF

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

The Python steps locate their folders from the repository structure and need no configuration.

## Running the Pipeline

Run the scripts in order from the `code/` folder:

```matlab
step0_setup
step1_discover
step2_preprocess
step3_erp
step4_timefreq
step5_rsa
```

Then from the shell:

```bash
python code/step6_decoding.py
python code/step7_specparam.py
```

Then back in MATLAB:

```matlab
step8_summary
```

## Citation

If using the ERP CORE dataset, cite the original ERP CORE paper:

Kappenman, E. S., Farrens, J. L., Zhang, W., Stewart, A. X., & Luck, S. J. (2021). *ERP CORE: An open resource for human event-related potential research*. NeuroImage, 225, 117465. https://doi.org/10.1016/j.neuroimage.2020.117465

If using the spectral parameterisation step, cite:

Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A. H., Wallis, J. D., Knight, R. T., Shestyuk, A., & Voytek, B. (2020). *Parameterizing neural power spectra into periodic and aperiodic components*. Nature Neuroscience, 23, 1655–1665. https://doi.org/10.1038/s41593-020-00744-x

## License

The code is released under the MIT License.
