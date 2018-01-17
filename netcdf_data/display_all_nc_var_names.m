load('netcdfstruct.mat')
var_nc_names=zeros(length(FINFO.Variables),100);

for n=1:length(FINFO.Variables)
     var_nc_names=FINFO.Variables(n).Name;
     disp(var_nc_names)
end