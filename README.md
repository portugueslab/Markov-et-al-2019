# Support code for Markov et al.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5147934.svg)](https://doi.org/10.5281/zenodo.5147934)


This repository contains the Python and MATLAB scripts that have been used to produce and analyze the data presented in Markov et al, 2021 (preprint [here](https://www.biorxiv.org/content/10.1101/2020.02.12.945956v1)).

The repo is organized in the following subfolders:
 - **experiment_control**: scripts to run all the behavioral and imaging experiments, based on [the Stytra library](https://github.com/portugueslab/stytra), and they provide a self-contained description of the behavioral paradigms described in the paper. To test them, you can refer to [the Stytra documentation](http://www.portugueslab.com/stytra/) for software installation and a description of the hardware. These protocols were run using the [0.8.5 Stytra version](https://github.com/portugueslab/stytra/commit/28171789a576f835021d4e0314ac923092aa3acc).
 - **preprocessing**: python notebooks and scripts to run all the preprocessing of the imaging data, using either the [fimpy library](https://github.com/portugueslab/fimpy) or the [suite2p library](https://suite2p.readthedocs.io/en/latest/).
 - **analysis**: MATLAB code for the analyses in the paper. Divided in three subsections:  
    1. **MATLAB_functions** (utility functions)
    2. **behavioral_analysis**
    3. **feedback_control_model**
    4. **imaging_analysis**


# Software
 - **MATLAB code**: all MATLAB code was developed using MATLAB 2018b
 - **Python code**: all Python scripts has been tested and run on Python 3.7. Stytra scripts were run using [Stytra 0.8.5](https://github.com/portugueslab/stytra/commit/28171789a576f835021d4e0314ac923092aa3acc). Suite2p calcium imaging data preprocessing was run using [suite2p 0.7.1](https://github.com/MouseLand/suite2p/releases/tag/0.7.1). Custom calcium imaging data preprocessing was done using code then moved to this [fimpy repository](https://github.com/portugueslab/fimpy).


# Installation and instructions to use
For replication of experiments using Stytra, you can refer to the [Stytra installation guide](https://www.portugueslab.com/stytra/userguide/0_install_guide.html) and the rest of Stytra documentation. 

For reproducing MATLAB-based analyis, you will need the data [uploaded with the manuscript](http://doi.org/10.5281/zenodo.5052785). 
Note that many functions use custom-written utility functions that can be found [here](https://github.com/portugueslab/markov-et-al/tree/master/analysis/MATLAB_functions), so please use `addpath(genpath(path/to/utility functions))` in order to reproduce the analysis.

### Behavior analysis
All behavioral data was acquired with Stytra and analyzed using custom-written MATLAB code [here](https://github.com/portugueslab/markov-et-al/tree/master/analysis/behavioral_analysis). This code uses raw compressed data uploaded with the manuscript and returns a master structure containing the processed information for every fish.
This information includes time points of swimming bouts, parameters of each bout (such as bout duration, peak tail beat amplitude, power profiles, duration of subsequent interbout, etc.), parameters for each trial (such as average bout/interbout duration, number of bouts, etc.) and metadata for each fish. 

### Imaging preprocessing
All imaging data was acquired and pre-processed using custom-written Python software or the published package [suite2p](https://github.com/MouseLand/suite2p/releases/tag/0.7.1). The pre-processing included alignment and detection of ROIs. We provide the code used for image preprocessing, but not the raw imaging data for size constraints. 

### Imaging analysis
All subsequent analysis was performed in MATLAB using custom-written code in [here](https://github.com/portugueslab/markov-et-al/tree/master/analysis/imaging_analysis). It is organized in three subfolders: 
 - whole_brain_imaging_inetgrators (experiment from Fig. 3)
 - PC_imaging_long_term_adaptation (experiment from Fig. 6)
 - whole_brain_imaging_long_term_adaptation (experiment from Fig. 7)
In each subfolder, there is a main filed called `[experiment_name]_analysis.mat`, and other files which are called by the main file. Each of the files performs a certain analysis step, such as anatomical registration of the ROIs to a common reference brain, computing triggered averages, sensory and motor scores, bar codes, time constants, etc., and saves the outcome of this step in a separate file. Note that we have only uploaded registered coordinates of ROIs, so reproduction of the registration step is not possible.

### Feedback control model
MATLAB code used for the feedback control model can be found at [here](https://github.com/portugueslab/markov-et-al/tree/master/analysis/feedback_control_model) and has a separate [readme.txt](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/readme.txt) file. 

### Figures generation
Finally, [Markov_et_al_2021_Make_Figures.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/Markov_et_al_2021_Make_Figures.m) takes analyzed data and produces all panels used in the Maniscript. 







For any question about the code, you can write at ruben.portugues@tum.de.
