function [wf] = compute_wf_v3 (t, taus, lat)

wf=t./(10.*(1-exp(-(lat-0.22)./taus)));

