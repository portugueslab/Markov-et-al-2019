function [wm] = compute_wm_v3 (taus, lat)

wm=1./(1-exp(-lat./taus));