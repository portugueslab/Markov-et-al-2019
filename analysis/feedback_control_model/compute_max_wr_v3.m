function [wr] = compute_max_wr_v3 (taus, lat)

wr=(-exp(-lat./taus))./(3.2.*(exp(-lat./taus)-1));