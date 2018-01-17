function [BS] = BS_add_filament(BS,Gamma,I,dGamma,ha)
%---------------------------------------------------
%  NAME:      BS_add_filament.m
%  WHAT:      Adds a filament to the BS analysis.
%  REQUIRED:  BS Toolbox 20150407
%  AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%  COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).
%
%  USE:
%  [BS] = BS_add_filament(BS,Gamma,I,dGamma)
%
%  INPUTS:
%    BS      = BS data structure
%    Gamma      = Filament points coordinates (x,y,z), one point per line [m,m,m]
%    I          = Filament current (flows from first point towards last point) [A]
%    dGamma     = Filament max discretization step [m]
%
%  OUTPUTS:
%    BS      = Updated BS data structure
%      BS.Nfilament              = Number of filaments
%      BS.filament(*).*          = Filament structure
%      BS.filament(*).Gamma      = Filament points coordinates (x,y,z), one point per line [m,m,m]
%      BS.filament(*).I          = Filament current (flows from first point towards last point) [A]
%      BS.filament(*).dGamma     = Filament max discretization step [m]
%----------------------------------------------------

n = BS.Nfilament+1;
BS.filament(n).Gamma = Gamma;
BS.filament(n).I = I;
BS.filament(n).dGamma = dGamma;
BS.Nfilament = n;

if nargin==5 && ~isempty(ha)
    %Plot P (where there is a current source)
    plot3(ha,Gamma(1,:),Gamma(2,:),Gamma(3,:),'.-r')
end
