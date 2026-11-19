# Standardizing 128-Channel Resting-State EEG Preprocessing and Spectral Parameterization

MATLAB pipeline for cleaning high-density (128-channel) resting-state EEG recordings and decomposing power spectral density (PSD) into broadband aperiodic (1/f) and narrowband periodic components using FOOOF / SpecParam.

## Repository Files

* **`NSCI270_Preprocessing.m`**: Automated EEG signal cleaning workflow (1 to 80 Hz FIR bandpass filtering, 58 to 62 Hz notch filtering, 200 Hz resampling, two-stage spherical spline bad-channel interpolation, common average re-referencing, and Extended Infomax ICA with ICLabel artifact rejection).
* **`NSCI270_FOOOF.m`**: Welch PSD estimation (2-second windows, 50% overlap), FOOOF / SpecParam model fitting (1 to 40 Hz), extraction of aperiodic exponents ($\chi$) and band-specific linear beta power (Low Beta: 13 to 16 Hz, Mid Beta: 16 to 20 Hz, High Beta: 20 to 30 Hz), and regional Spearman correlation analyses.

## Requirements

* **MATLAB** (tested on macOS)
* **EEGLAB** (with **ICLabel** plugin and EGI `.mff` import support)
* **FOOOF / SpecParam** (MATLAB wrapper or Python `specparam` environment)

## Usage

1. Update the input and output file paths at the top of `NSCI270_Preprocessing.m` and run it on raw resting-state recordings.
2. Point `NSCI270_FOOOF.m` to the cleaned EEG files to compute Welch PSDs, fit aperiodic and periodic models, and export summary tables and figures.# eeg_preprocessing_fooof_pipeline_nsci270
