# Feedback control model

The core of the model is [brain_iteration_v3.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/brain_iteration_v3.m)

It computes two variables:
1. brain_state - current state of the model nodes: forward and reverse velocity sensors, velocity and motor integrators, and motor output generator.
2. swim - binary swimming variable

... based on:
1. Previous brain_state
2. Previous swim
3. Current grating speed (mm/s)
4. Set of 8 model parameters

Another useful function is [exp_iteration_v3.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/exp_iteration_v3.m)
It computes current grating speed depending on behavior of the model and reafference condition. 
Reafference condition is defined as a 1 x 5 array with the following columns:
1. gain
2. lag
3. shunted (true) or non-shunted lag (false)
4. gain drop start (measured in frames with respect to the bout onset)
5. gain drop end

These two functions are combined in an experiment program that defines an artificial experimental protocol.
Protocols used in this study can be found in the model_protocols folder:
1. [model_one_bout_trial_v3](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_protocols/model_one_bout_trial_v3.m) - grating starts moving 0.7 s after beginnig of the protocol
and stops 0.4 s after bout offset. The trial stops 1.1 s after bout offset. This protocol was used to present small pictograms of model nodes in Fig. 2b.
2. [model_v3_real_trial.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_protocols/model_v3_real_trial.m) - same as trials used in the acute reaction experiment. It was used to show model traces in an example trial (Fig. S1).
3. [model_short_trial_v3.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_protocols/model_short_trial_v3.m) - grating starts moving 0.3 s after beginnig of the protocol. The trial terminates at the onset of the 2nd bout.
This protocol was used for fitting the model to real data (Fig. 2c).

Finally, [model_compute_parameters_v3](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_compute_parameters_v3.m) simply computes second bout and interbout duration after [model_short_trial_v3.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_protocols/model_short_trial_v3.m) (Fig. 2c).

[model_v3_individual_fish_genetic_fitting.m](https://github.com/portugueslab/markov-et-al/blob/master/analysis/feedback_control_model/model_v3_individual_fish_genetic_fitting.m) was used to fit model parameters to all individual fish tested in the acute reaction experiment.
