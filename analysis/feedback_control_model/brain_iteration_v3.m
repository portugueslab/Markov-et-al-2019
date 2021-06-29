function [swim, brain_state] = brain_iteration_v3(swim, brain_state, grspeed, par)
%
% Inputs:
%
% 1. swim - a binary variable that tells if fish swam before this iteration
%
% 2. brain_state - previous state of the brain
%       brain_state(1) - activity of forward motion sensor
%       brain_state(2) - activity of reverse motion sensor
%       brain_state(3) - activity of sensory integrator
%       brain_state(4) - activity of motor output generator
%       brain_state(5) - activity of motor integrator
%
% 3. grspeed - current grating speed
%
% 4. par - parameters of the model
%       par(1) - wf      - weight between forward motion sensor and sensory integrator
%       par(2) - wr      - weight between reverse motion sensor and sensory integrator
%       par(3) - dt/taus - time constant of sensory integrator
%       par(4) - wi      - weight between motor integrator and motor output generator
%       par(5) - ws      - weight of feed-forward self-excitation of motor output command cell
%       par(6) - t       - threshold of motor output command
%       par(7) - wm      - weight between motor output command cell and motor integrator
%       par(8) - dt/taum - time constant of motor integrator
%
%
% Outputs:
%
% 1. swim - a binary variable that tells if fish swims after this iteration
% 
% 2. brain_state - state of the brain after this iteration

% forward motion sensor (positively rectified grating speed)
brain_state(1)=max(grspeed,0);

% reverse motion sensor (negatively rectified grating speed)
brain_state(2)=-min(grspeed,0);

% sensory integrator (leaky integrator with saturation at 1)
brain_state(3)=max(min(par(3)*par(1)*brain_state(1)-par(3)*par(2)*brain_state(2)-(par(3)-1)*brain_state(3),1),0);

% motor output generator (activated by sensory integrator and inhibited by motor integrator)
brain_state(4)=max(brain_state(3)-par(4)*brain_state(5),0);

% swim (fish swims if motor output generator + self-excitation is greater than swimming threshold)
swim=brain_state(4)+par(5)*swim>par(6);

% motor integrator (leaky integrator with saturation at 1)
brain_state(5)=min(par(8)*par(7)*swim-(par(8)-1)*brain_state(5),1);