# Krok-etal_2023_DA-ACh

The functions and script in this repository were developed to process, analyzem and visualize data acquired during fiber photometry and electrophysiological experiments presented in Krok et al. 2023.

This pipeline includes the following processing / analysis:

* Fiber Photometry
  * Filtering, Baselining, Demodulation
* Movement
  * Velocity Calculation, Movement Onset and Offset Alignment
* Reward
  * Reward Delivery and Lick Onset Alignment
* Figures
  * Plotting figures in Krok et al. 2023


## Dependencies:

The pipeline was written using MATLAB2018a

Tritsch Lab Toolbox available at: https://github.com/TritschLab/TLab_Toolbox

Tritsch Lab Photometry Pipeline available at: https://github.com/pratikmistry96/Photometry-Pipeline

Buzcode available at: https://github.com/buzsakilab/buzcode

Chronux Toolbox available at: http://chronux.org

Please download these toolboxes and pipelines and add to your path.

## Installation:

Download the repositories and add them to your MATLAB path

## Steps:

1. Convert h5 files to .mat files

      From the TLab Toolbox, run **convertH5** to convert h5 files from wavesurfer into the format required for analysis, .mat

       >> convertH5

2. Edit Parameters
      
      From the Photometry Pipeline, open the **processParams** file and adjust the parameters accordingly

       >> edit processParams

      Save the file with a new name

3. Process Data

      From the Photometry Pipeline, run **processData**. Select the desired parameter file (from step 2) for the experiment. The processed data file will overwrite the .mat file generated in Step 1.

       >> processData
    
4. Extract Data into single structure

      The **extractBeh** script will extract processed signals from multiple .mat files.

       >> extractBeh


## Author:

Name: Anne Krok

Email: ack466@nyu.edu
