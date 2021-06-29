function [par_ranges] = compute_par_ranges_v3 (max_taus, max_taum, max_wi, lat_ranges, tlat, tmot, par_ranges)
% computes the most conservative parameter ranges 
% see labfolder entry as of 08.05.2019 and 09.05.2019 in Daniil LabBook for details

% function inputs:
% max_taus - maximal time constant of the sensory integrator [s]
% max_taum - maximal time constant of the motor integrator [s]
% max_wi - maximal weight between motor integrator and motor output generator
% lat_ranges - range of allowed latency to initiate swimming after grating onset [s]
% tlat = 0.01 - minimal allowed time for comlete discharge of the sensory integrator (from 1 to 0) during a bout at gain 0.66
% tmot = 0.01 - minimal allowed time of complete saturation of the motor integrator (from 0 to 1) during a bout [s]
% par_ranges - optional ranges (useful for the fitting algorythm). If this is defined, the function will only set dependent ranges

% output is n x 2 array, where n - is number of model parameters, 
% 1st column is the min limit, 2nd column is the max limit
% par(1) - wf      - weight between forward motion sensor and sensory integrator
% par(2) - wr      - weight between reverse motion sensor and sensory integrator
% par(3) - taus    - time constant of sensory integrator [s]
% par(4) - wi      - weight between motor integrator and motor output generator
% par(5) - ws      - weight of feed-forward self-excitation of motor output command cell
% par(6) - t       - threshold of motor output command
% par(7) - wm      - weight between motor output command cell and motor integrator
% par(8) - taum    - time constant of motor integrator [s]

% create an array and write down independent ranges (if it wasn't defined)
if nargin == 6
    par_ranges=[...
        nan nan;...     % wf
        0   nan;...     % wr
        0   max_taus;...% taus
        0   max_wi;...  % wi
        0   nan;...     % ws
        0   1;...       % t
        0   nan;...     % wm
        0   max_taum;...% taum
        ];
end

% define ranges for wf
par_ranges(1,1)=compute_wf_v3(par_ranges(6,2), par_ranges(3,2), lat_ranges(2));
par_ranges(1,2)=compute_wf_v3(par_ranges(6,2), par_ranges(3,2), lat_ranges(1));

% define upper limit for wr
par_ranges(2,2)=compute_max_wr_v3(par_ranges(3,2), tlat);

% define  upper limit for ws
par_ranges(5,2)=par_ranges(6,2);

% define upper limit for wm
par_ranges(7,2)=compute_wm_v3(par_ranges(8,2), tmot);