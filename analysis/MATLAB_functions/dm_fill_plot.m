function [] = dm_fill_plot(time,mean_trace,ste_trace,col,alpha)
fill([time flip(time)],[mean_trace-ste_trace flip(mean_trace+ste_trace)],col,'edgecolor','none','facealpha',alpha);
plot(time,mean_trace,'color',col,'linewidth',1.5);
