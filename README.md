# README #

This document outlines all scripts used to process data for RNA Polymerase II is a Polar Roadblock to a Progressing DNA Fork

## MATLAB Version ##

This code was developed using MATLAB.

MATLAB Version: 24.2.0.2863752 (R2024b) Update 5

Operating System: Microsoft Windows 11

MATLAB R2024b can be downloaded from MathWorks (https://www.mathworks.com/products/matlab.html) and installed with license. MATLAB can be installed within 30 min.

## Examples ## 

Each script comes with its own respective demo dataset file. The script and its respective demo data file are required to be in the same folder to run. Scripts are written, such that the corresponding sample data is imported upon running.

## Instructions for Use ##

1. Download scripts from Github Repository using `Download ZIP` (approximate download time <30s)
2. Prepare the environment: 
- Place all required data files and scripts in the same directory
- Ensure all custom functions are accessible in the MATLAB path. 
3. Launch MATLAB and navigate to the directory containing the script
4. Run the script and follow prompts

## Scripts ## 
### Figure1and2_AnalysisPlotting.m ###

This MATLAB script plots and analyzes experimental single-molecule data containing DNA stretching and unzipping. This script allows users to measure the transcribed distance, the disruption force, the sliding distance, and plot the extension shift due to RNA-DNA hybridization for visualization.  

Time to run: The first run of this script can take a little longer than the average run, due to a parallel processing environment in calculate_velocity_position.m. Subsequent runs after initial parallel processing environment is launched will take less than a minute. Not counting user interaction time.

Files Required:
- HeadOn_Template_Theory.txt
- Figure 1 and 2_Sample Data_Trace_02366.dat
- calculate_velocity_position.m

User Interaction: For traces with large force deviation from DNA baseline
During execution, a message box will prompt you to zoom into a region of interest. Use the zoom tool to zoom into region of interest in the Figure 1 window. Once zoomed, press Enter to continue. You will then need to click on a point within the zoomed region to mark the end of sliding. 

ROI: The last region of populated data points prior to the return to the DNA unzipping baseline. 

Output Metrics: Displayed in the respective figure windows
- Transcribed distance – Figure 1 
- Displacement force – Figure 1 
- Sliding distance – Figure 5 

### Figure3and4_AnalysisPlotting.m ###

This MATLAB script plots and analyzes experimental single-molecule data containing DNA stretching, DNA unzipping, RNA Pol II backtracking, and DNA reannealing. This script allows users to measure the backtrack distance, and the ability/inability for DNA to reanneal.

Time to run: Less than a minute. Not counting user interaction time. 

Files Required:
- HeadOn_Template_Theory.txt
- Figure 3 and 4_Sample Data_Trace_02531.dat

User Interaction: 
During execution, a message box will prompt you to zoom into a region of interest. Use the zoom tool to zoom into region of interest in the Figure 4 window. Once zoomed, press Enter to continue. You will then need to click on a point within the zoomed region to mark the end of sliding. 

ROI: Minimum force prior to the return to the theoretical DNA baseline.

Output Metrics: Displayed in the respective figure windows
- Transcribed distance – Figure 1 
- Backtracked Distance - Figure 1 and 3
- Force Dip during Rezipping - Figure 4 
- Force at Ext = 2750nm during 2nd Unzip - Figure 2 

### Figure5_AnalysisPlotting.m ###  

This MATLAB script plots and analyzes experimental single-molecule data containing DNA stretching, DNA unzipping, RNA Pol II backtracking, and DNAP replication. This script allows users to measure the backtrack distance, the ability/inability for DNA to reanneal, and DNAP activity.

Time to run: Less than a minute. 

Files Required:
- HeadOn_Template_Theory.txt
- Figure 5_Sample Data_Trace_01502.dat

User Interaction: This script runs automatically. No manual input is required during execution. 

Output:
- Transcribed distance – Figure 1 
- Backtracked Distance - Figure 1 
- ssDNA available to be replicated - Figure 1 
- Number of base pairs replicated - Figure 2 


### calculate_velocity_position.m ###
Function that is used to calculate the velocity and smoothed position from the position vs time data. Used in (`Figure1and2_AnalysisPlotting.m`)
