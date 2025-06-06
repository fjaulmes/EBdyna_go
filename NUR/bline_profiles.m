function [ne_b,Te_b,Zeff_b]=bline_profiles(R,Z,ne,Te,Zeff,Rg,Zg)


ne_b=interp2(Rg,Zg,ne,R,Z,'*linear');
Te_b=interp2(Rg,Zg,Te,R,Z,'*linear');
Zeff_b=interp2(Rg,Zg,Zeff,R,Z,'*linear');


end