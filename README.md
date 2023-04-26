README for Supplementary material of Sabine C. Fischer, Simon Schardt, Joaquín Lilao-Garzón and Silvia Muñoz Descalzo  "The salt-and-pepper pattern in mouse blastocysts is compatible with signalling beyond the nearest neighbours."

All code is provided as Mathematica .nb. If a Mathematica licence is not available, the files can be used with the Wolfram Player, available at https://www.wolfram.com/player/.

"dataPreprocessing" contains individual folders for all experimental data sets analysed and the simulation data from the rule-based as well as ODE-based signalling models.

Each folder for the experimental data contains at least three subfolders:
- "DataFromMINS" contains the output from the image analysis software MINS (Lou et al. Stem Cell Rep. 2014)
- "EmbryoFeatures" contains the processed data in Mathematica's .mx format
- "Results" contains the output from the preprocessing steps. In particular, the files *Comp.mx that contain all the information for the neighbourhood analyses.

For the simulation results, the folder "simulations" contains 
- "DataFromPython" with the results from the simulations in Python of the signalling models (Schardt & Fischer, arXiv, 2022)
- "TissueFeatures" with the processed data in Mathematica's .mx format
- "Results" with the output from the preprocessing steps. In particular, the files *Comp.mx contain all the information for the neighbourhood analyses.

The folder "analysisCode" contains all calculations that use the complete experimental and the complete simulations data set. Notebooks "Neighbour*.nb" calculate mainly the ICM composition scatter plots from the "*Comp.mx" files. "MoranIndex*.nb" contain mainly the calculations of the Moran's index.

All further information regarding the analysis is contained in the Mathematica notebooks. (Results of longer calculations are saved as .mx files. The respective code lines are commented out and the resulting .mx files provided.)

All figures that were generated with Mathematica can be found in the subfolders "Figures". The remaining figures in the paper were generated with Excel based on the .xls files provided in the "Results" folders.

"dataJson" contains the processed data and the neighbour comparison files for data sets I-VI in .json format for compatibility with other analysis software. 


For questions and comments please contact sabine.fischer@uni-wuerzburg.de

The latest release is archived at Zenodo: 



