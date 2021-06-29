function regressor = trace2regressor(trace,my_kernel,time_be,time_im)
regressor=conv(trace,my_kernel);
regressor(length(trace)+1:end)=[];
regressor=interp1(time_be,regressor,time_im);
% if any(regressor>0)
%     regressor=regressor/max(regressor);
% end
temp2=find(isnan(regressor));
regressor(temp2)=[];
regressor=zscore(regressor);
for i2=1:length(temp2)
    regressor=[regressor(1:temp2(i2)-1) nan regressor(temp2(i2):end)];
end
regressor(temp2)=nan;