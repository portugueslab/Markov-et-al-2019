function [col] = dm_lab2rgb(col)
col=min(max(lab2rgb(col),0),1);