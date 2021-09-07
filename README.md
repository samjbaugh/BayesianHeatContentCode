# README


Code and data for "Hierarchical Bayesian Modeling of Ocean Heat Content and Its Uncertainty" which is under review at AOAS as of 9/7/21. Included is the source code and documentation for the BayeysianOHC package as well as scripts for reproducing the results, tables, and figures in the paper. Among other functionality this includes computing the initial configuration for the sampler, sampling the posterior distribution using  MCMC for either the full global model or for sub-regions, post-processing the posterior samples, fitting comparison models including replication of the Levitus et al. (2012) approach,  computing cross-validation scores, and evaluating the cylindrical distance convolutions. 

In addition, the Argo data, initial configuration output, Vecchia approximation structures, 20,000 posterior samples for the main model fit, and data created in the course of the analysis are included. 

This repository contains the following directories:

* AnalysisScripts: Scripts for conducting the analysis. See the included README for documentation of each script's inputs and outputs. Each script is self-contained but may require functions from the tidyverse package.
* BayesianOHC: Source code for BayesianOHC pacakge. See BayesianOHC_0.1.0.pdf for documentation.
* FigureScripts: Scripts for reproducing figures in the paper.
* MCMC_Input: Processed data needed for initializing the MCMC sampler. Includes:
	* argo_data_january.RData
	* argo_data_subset.RData
	* grouping_list.RData
	* initial_parameters.RData
	* knot_points.RData
* MCMC_Output: Output of the MCMC sampler for the fully non-stationary run. Includes the "Output_20k_run" directory with the following contents: 
	* map_configuration.RData
	* mu_slope_df.RData
	* ohc_df.RData
	* posterior_intervals.RData
	* sample_matrices.RData
	* trend_resamp_intonly.RData
	* trend_resampled.RData*
* RawData: Unprocessed data; includes the following:
	* mask.RData
	* argo_object_01.RData*
* ValidationData: Data used for the cross validation study in the paper. See the README in the AnalysisScripts folder to see where each file is generated and used. lowodfs and map_params (which are large since they contain veccmat objects) are not tracked in order to keep the size down.
	
(*) Files that are in .gitignore in order to keep the repository from being too large
