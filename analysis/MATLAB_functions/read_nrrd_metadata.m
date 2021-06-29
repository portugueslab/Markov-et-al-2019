function [sz, rez] = read_nrrd_metadata (str_nrrd)
[~, rez]=nrrdread(str_nrrd);
sz=strread(rez.sizes);
sz=[sz(2) sz(1) sz(3)];
rez=rez.spacedirections;
brackets=find(rez=='(');
commas=find(rez==',');
rez=[str2double(rez(brackets(1)+1:commas(1)-1)) str2double(rez(commas(3)+1:commas(4)-1)) str2double(rez(commas(6)+1:end-1))];
rez=[rez(2) rez(1) rez(3)];