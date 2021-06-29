function [ca_tau] = genotype2tau(genotype)
if contains(genotype,'GCaMP6s')
    ca_tau=1.8;
elseif contains(genotype,'GCaMP6f')
    ca_tau=0.4;
elseif contains(genotype,'GCaMP6sef05') || contains(genotype,'GCaMP6fef05')
    ca_tau=0.65;
else
    error ('Unknown genotype')
end