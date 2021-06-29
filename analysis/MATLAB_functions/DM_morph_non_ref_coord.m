function [ref_coord]=DM_morph_non_ref_coord(non_ref_coord,xfrom_dir_path,pathname)
if pathname(end)~='\'
    pathname=[pathname '\'];
end
vlist(:,1)=non_ref_coord(:,2);
vlist(:,2)=non_ref_coord(:,1);
vlist(:,3)=non_ref_coord(:,3);
f_list='temp.list';
dlmwrite(fullfile(pathname,f_list),vlist,'delimiter','\t','newline','pc');
str_list_in=strcat('<',pathname,f_list,'>');
str_list_out=strcat(pathname,strrep(f_list,'.list','_out.list'));
dos(['C:\cmtk_files\bin\streamxform' ' ' '-- --inverse' ' ' xfrom_dir_path ' ' str_list_in ' ' str_list_out]);
fid=fopen(fullfile(pathname,strrep(f_list,'.list','_out.list')));
myline=fgetl(fid);
rr=0;
while(ischar(myline))
    rr=rr+1;
    substrings=regexp(myline,' ','split');
    ref_coord(rr,2)=str2num(substrings{1});
    ref_coord(rr,1)=str2num(substrings{2});
    ref_coord(rr,3)=str2num(substrings{3});
    if(numel(substrings)==4)
        ref_coord(rr,4)=0;
    else
        ref_coord(rr,4)=1;
    end
    myline=fgetl(fid);
end
fclose(fid);
ref_coord(:,4)=[];