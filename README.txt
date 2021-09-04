Author: Jamie Cooper Hockey
Date: 4/9/2021

Directory contents:

DNNGP_CODE:
	- data:
		- SYNTH_DATA: synthetic data used for dnngp
		- TEMP_DATA: temperature data used for dnngp
	- dnngp:
		- Makefile
		- nngp.cpp
		- pfiles-realdata: input files for 3 different chains for temperature data
		- pfiles-syntheticdata: input files for 3 different chains for synthetic data
	- libs: library files used for nngp.cpp
NETemp:
	- dynLM_tempdata: code for running dynamic linear model on temperature data
	- NEtemp_data_preprocessing: pre-processing of temperature data
	- processing_results_dnngp_realdata: post MCMC processing of dnngp temperature data outcomes
	- prediction_tempdata_dnngp: code for posterior prediction from dnngp
synthetic_data:
	- dynLM_synthetic: code for running dynamic linear model on synthetic data
	- generate_synthetic_data: synthetic data generation code
	- processing_results_dnngp_synth: post MCMC processing of dnngp synthetic data outcomes
	- prediction_synthetic_dnngp: code for posterior prediction from dnngp
Visualizations:
	- covariance_matrix_visualization_data: code to generate data used for covariance matrix visualization
	- covariance_matrix_approximation: code to compute and visulaize inverted covariance and approximate covariance matrices
	- data_visualization: code to visualize temperature data


Overview:
dnngp code slightly adapted from "response-matern" code taken from (Finley et al., 2019). dnngp code was run on the "Hamilton" Durham University HPC.



Citations:

Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., and Banerjee, S.
(2019), “Efficient Algorithms for Bayesian Nearest Neighbor Gaussian Processes,” Journal
of Computational and Graphical Statistics, 28, 401–414, pMID: 31543693.