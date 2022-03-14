# Sensitivity-Analysis-of-a-CSG-climate-model

For running any of the scrips make sure that this and the folder "2 Raw Data" are added to the Matlab Path. In this folder you can find the code used for the calculation of sensitivity indices; the content is divided as:

	3-0 CSG Model. Developed by Zhou (2020). 

	3-1 One-at-a-time SA. Containing the code for geneerating the normalized sensitivity indices. 
	
	3-2 Multiple linear regression SA. Here are two scripts, one that calculates the standardise regression coefficients (SRC) using two methods:
	
		linearRegMatrices.m By matrix algebra
		scatterAndLinearReg.m By fitlm Matlab function. It also plots the corresponding scatterplots (Figure 4.5). 
	3-3 Sobol'indices SA. This folder contains use functions to generate samples using MonteCarlo and Latin Hypercube, and two Matlab scripts:
		generateSobolMatrices.m to evaluate the model for matrices A, B, C and D and generate their corresponding outputs. 
		calculateSobolIndices.m which calculates the first, second and total order Sobol' indices, and generates Figure 4.8.
	3-4 Plotting. Here you can find the scripts to generate the figures of this thesis. 
		usefulPlots.m generates Figure 1.1
		prettyPlots.m generates Figures 4.1, 4.2, 4.3, 4.4, 4.6, 4.7, and 5.1 
