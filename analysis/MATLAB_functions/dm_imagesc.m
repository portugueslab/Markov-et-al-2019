function dm_imagesc(time,data)
dt=time(2)-time(1);
num_traces=size(data,1);
set(gca,'xlim',[time(1)-dt time(end)]+dt/2,'ylim',[0 num_traces]+0.5,'ydir','reverse','layer','top');
imagesc([time(1) time(end)],[1 num_traces],data,'AlphaData',double(~isnan(data)));
line([0 0],[0 num_traces]+0.5,'color','k','linestyle',':');
c1=[40 75 62];
c2=[90 0 0];
c3=[40 16 -69];
l = [linspace(c3(1), c2(1), 32) linspace(c2(1), c1(1), 32)];
a = [linspace(c3(2), c2(2), 32) linspace(c2(2), c1(2), 32)];
b = [linspace(c3(3), c2(3), 32) linspace(c2(3), c1(3), 32)];
my_colormap=[l;a;b]';
my_colormap(32,:)=[];
my_colormap=dm_lab2rgb(my_colormap);
c_lim=get(gca,'clim');
c_lim=max(abs(c_lim));
c_lim=[-c_lim c_lim]*0.9;
set(gca,'clim',c_lim);
colormap(my_colormap);
    