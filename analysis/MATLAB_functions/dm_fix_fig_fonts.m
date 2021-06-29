function [] = dm_fix_fig_fonts(h,font_size)
drawnow
if nargin==0
    h=gcf;
end
if nargin<2
    font_size=10;
end
set(h,'InvertHardcopy','off','color',[1 1 1]);
all_axes=findall(h,'type','axes');
set(all_axes,'fontsize',font_size,'fontweight','bold','FontName','Arial','layer','top','tickdir','out');
for i=1:length(all_axes)
    this_xcol=all_axes(i).XColor;
    if isequal(this_xcol, [0.15 0.15 0.15])
        all_axes(i).XColor=[0 0 0];
    end
    this_ycol=all_axes(i).YColor;
    if isequal(this_ycol, [0.15 0.15 0.15])
        all_axes(i).YColor=[0 0 0];
    end
end
all_polaraxes=findall(h,'type','polaraxes');
set(all_polaraxes,'fontsize',font_size,'fontweight','bold','FontName','Helvetica','layer','top');
all_colorbars=findall(gcf,'type','colorbar');
set(all_colorbars,'fontsize',font_size,'fontweight','bold','FontName','Helvetica','tickdir','out');
