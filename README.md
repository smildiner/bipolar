# Project Name: Bipolar

This repository contains R code to reproduce the results of "Mildiner Moraga, S., Bos, F. M., Doornbos, B., Bruggeman, R., van der Krieke, L., Snippe, E., & Aarts, E. (2023, January 31). Evidence for severe mood instability in patients with bipolar disorder: Applying multilevel hidden Markov modelling to intensive longitudinal ecological momentary assessment data. <https://doi.org/10.31234/osf.io/egp82>".

## The repository includes three main folders:

 - `data`: empirical data folder (not included due to privacy issues).
 - `outputs`: contains the folders `figures` (figures included in the manuscript) and `res` (each of the outputs for the simulated data).
 - `R`: R scripts to fit Bayesian multilevel HMM on empirical data, analyse main results, run simulation, decode simulation results, and analyse simulation results, and necessary utility functions.
 - `results`: main analysis results from fitting the Bayesian multilevel HMM on the empirical data; contains the main object generated as output after fitting the multilevel HMM.

## Reproducing Results

To reproduce the results of the study, follow these steps:

1. Download the code in the `R` folder and the empirical data (upon request from researchers).
2. Run the R scripts in the following order:
  1. `1_fit_mhmm.R`: fit the Bayesian multilevel HMM on the empirical data.
  2. `2_main_results.R`: main results from fitting the Bayesian multilevel HMM on the empirical data.
  3. `3_supp_results.R`: supplementary results from fitting the Bayesian multilevel HMM on the empirical data.
  4. `4_run_simulation.R`: run the simulation.
  5. `5_analyse_simulation_results.R`: decode and analyse the simulation results.

Notice that the simulation is computationally intensive; it was run in the [Dutch National Supercomputer Snellius](https://www.surf.nl/en/dutch-national-supercomputer-snellius) using a single node with 128 cores.

## Contact

The empirical data cannot be found in this repository due to privacy issues. However, the data can be made available upon reasonable request from researchers. Please contact <f.m.bos01@umcg.nl> for more information.

For any further questions, please contact <s.mildinermoraga@uu.nl>.