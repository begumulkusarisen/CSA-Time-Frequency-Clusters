# K-Clustering-for-EEG-Time-Frequency-Coherence-data-from-Brainstorm
The repository contains MATLAB code for the k-clustering of the EEG functional dynamic connectivity matrices derived from Brainstorm. The script was made with the help of Gemini 3 Flash.

This script processes connectivity matrices (e.g., coherence data) across multiple subjects to identify recurring brain states using K-means clustering. 
As our Brainstorm analysis was based on the Schafer Atlas, this script further focuses on specific functional networks (Default, Control, and Dorsal Attention), flattens the connectivity features, and calculates temporal metrics like Transition Probability Matrices (TPM) and Dwell Times.

Workflow:

# Data Extraction & Preprocessing 
File Loading: Automatically searches for connectivity files with the pattern *connectn_cohere_*.ROI Mapping: Parses ROI labels to identify which nodes belong to the Default, Cont, and DorsAttn networks based on naming conventions.3D Reconstruction: Converts flattened connectivity vectors back into square symmetric matrices for each time point. Feature Flattening: Extracts the upper triangle of the matrices (to avoid redundant edges) and calculates absolute values for feature input.

# Group-Level Clustering (K-Means)
The script stacks data from all subjects into a massive group matrix for each network.
K-Means Clustering (k=5): Groups time points into 5 distinct "states" based on connectivity patterns using squared Euclidean distance.Subject Mapping: Re-segments the group cluster labels back to individual subjects to preserve temporal order.

# Temporal Dynamics Analysis 
Transition Probability Matrices (TPM): Calculates the probability of switching from State $A$ to State $B$.Dwell Time: Computes the average consecutive time (number of frames) a subject remains in a specific state before transitioning.

# Visualization & ExportHeatmaps 
Generates group-averaged TPM heatmaps showing the most common state transitions. Excel Export: Saves a multi-sheet Excel file (Subject_Dwell_Times_Results.xlsx) containing individual dwell times for every subject across all networks.

# Requirements: 
1) MATLAB (Statistics and Machine Learning Toolbox required for k-means).
2) .mat files containing a TF matrix (connectivity) where the Atlas or RowNames variable exists for ROI identification.

