function [maps]=split_theta_map(maps)
%%split_theta_map Creates second theta map in maps struct.
% One needs a secondary theta map if one wants to find theta by
% interpolation. This is due to interpolation error likely for points
% with angle higher than max(theta), but lower than 2*pi (e.g. 2*pi-eps).
% Interpolation in nodes with 2*pi difference is avoided by making sure the
% transition is smooth (making use of negative numbers).


%% Check if condition if maps-struct
if isfield(maps,'theta_normal_XZ') || isfield(maps,'theta_phase_shift_XZ')
    warning('Theta maps already split')
    return
elseif ~isfield(maps,'theta_XZ')
    error('Theta map not found improperly split')
end

%% Add the normal and phase shifted map
maps.theta_normal_XZ=maps.theta_XZ;
maps.theta_phase_shift_XZ=mod(maps.theta_XZ+3.14,2*pi)-3.14; % Add pi, mod and substract pi to overcome zero crossing problem

%% Remove the obsolete 'normal' map
maps=rmfield(maps,'theta_XZ');

end