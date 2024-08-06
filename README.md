# SASKit-ML Pipeline Description
## Introduction
This document explains the SASKit-ML Pipeline, the machine learning component of the published SASKit protocol detailed in Henze et al. 2020 #REF and tested in Palmer et al. 2024. This code can be used to replicate the main results from Palmer et al. 2024 for the HMP-T2D, HMP-IBD and RA-MAP data sets. The same pipeline will also work for user downloaded TCGA data, although these are not provided in the repository due to size constraints. For help in reproducing the TCGA results please contact the authors.
## Pipeline
The SASKit-ML pipeline draws on three main data types, clinical blood markers, transcriptomics and proteomics, in addition to age and sex. These data sources are used to calculate up to 25 integrated features (plus age and sex) which are then used in the training of a random survival forest model predicting health deterioration (according to disease specific endpoint definitions).
Clinical markers are used directly as features without modification. Both transcriptomic and proteomic data are integrated to identify the most active pathways and subnetworks in the patients compared to healthy controls, under the assumption that expression changes in active subnetworks and pathways associated with the disease of interest relative to healthy individuals will be predictive of more severe disease and worse long-term outcomes.

The omics data are integrated in two ways to capture different aspects of the expression changes in the disease of interest. 

The first method is by using the popular GSVA tool to calculate “pathway activation scores” for KEGG pathways. This is achieved by calculating the GSVA score for each pathway for each patient/control and then performing a “differential activation analysis” using limma to identify those most strongly differentially activated pathways in the patients relative to controls. The top-5 of these pathways are then chosen, and the activation scores calculated for individual patients of these 5 pathways then used as the GSVA features in the random survival forest model.

The second method is by using ExprEssence, a method for identifying active subnetworks for a given network type (in this case, we use STRING physical protein-protein interactions), based on identifying regions in the network where multiple genes/proteins display large expression changes in patients compared to controls. This is achieved by calculating “linkscores” for all the edges in the STRING graph, then restricting the graph to only the top-50 edges ranked by this linkscore, thus reducing the larger graph into small subnetworks. The full linkscore calculation method is described in the ExprEssence paper (#ref). Once these subnetworks have been determined, the average linkscore for each is calculated and these are used as features in the random survival forest model. In the event that more than 5 subnetworks are identified, the subnetworks are ranked by their absolute average linkscore, and the average linkscores of only the top-5 by this ranking are used as features in the random survival forest model.

In summary, the total features for the random survival forest are:

* Up to 5 clinical markers (neutrophil–lymphocyte ratio, fibrinogen, high-sensitive C reactive protein, albumin and PAI-1).
* Up to 5 transcriptomic GSVA features.
* Up to 5 transcriptomic ExprEssence features.

If high throughput proteomics are available:

* Up to 5 proteomic GSVA features.
* Up to 5 proteomic ExprEssence features.

These features are then used as features in a 10-fold cross-validated random survival forest model predicting the time to a disease specific health deterioration outcome.
## Replicating the Manuscript Results
To run the analysis from the paper, simply download the folder "Analysis" and run "saskitMLScript.R" from within it.
