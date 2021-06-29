# Markov et al 
This repository contains the Python and MATLAB scripts that have been used to produce and analyze the data presented in Markov et al, 2021 (preprint [here](https://www.biorxiv.org/content/10.1101/2020.02.12.945956v1)).

The repo is organized in the following subfolders:
 - **experiment_control**: scripts to run all the behavioral and imaging experiments, based on [the Stytra library](https://github.com/portugueslab/stytra), and they provide a self-contained description of the behavioral paradigms described in the paper. To test them, you can refer to [the Stytra documentation](http://www.portugueslab.com/stytra/) for software installation and a description of the hardware. These protocols were run using the [0.8.5 Stytra version](https://github.com/portugueslab/stytra/commit/28171789a576f835021d4e0314ac923092aa3acc).
 - **preprocessing**: python notebooks and scripts to run all the preprocessing of the imaging data, using either the [fimpy library](https://github.com/portugueslab/fimpy) or the [suite2p library](https://suite2p.readthedocs.io/en/latest/).
 - **analysis**: MATLAB code for the analyses in the paper. Divided in three subsections:  
    1. **MATLAB_functions** (utility functions)
    2. **behavioral_analysis**
    3. **feedback_control_model**
    4. **imaging_analysis**


For any question about the code, you can write at ruben.portugues@tum.de.
