function [ pphi_kin ] = get_pphi_kin(input,x,v_plus_1,psi)
%get_pphi_kin Determines pphi_kin/eV for 2D and 3D maps. Note the full maps
%need to be loaded in `maps'.
%   Detailed explanation goes here

global par const

%% Pphi calculation
if par.APPLY_3D
    [field_3D] = find_3D_Afield(x,{'Aphi'});
    
    pphi_kin=...
        (input.m/const.eV)*x(:,1,:).*v_plus_1(:,3,:)...
        +input.Z*	(-psi...
        +x(:,1,:).*field_3D.Aphi);
else
    %% 2D fields
    pphi_kin=...
        (input.m/const.eV)*(x(:,1,:)).*v_plus_1(:,3,:)...
        -input.Z*psi;
end
end

