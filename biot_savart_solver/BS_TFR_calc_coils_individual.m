function [TF_map] = BS_TFR_calc_coils_individual(par,TF_map,TF_coils,Z_correction)
%BS_RMP_calc_coils_individual Summary of this function goes here
%   Detailed explanation goes here

%% Calculate B-field for each coil:
BS(1).Nfilament=1;

disp('************************************************')
disp('Starting standard B-field calculation of each coil:')
disp('************************************************')
disp(['length(TF_coils)= ' num2str(length(TF_coils))])

%% Correction for offset magnetic axis
TF_map.Z=TF_map.Z+Z_correction;   % Reduce Z express the vertical distance from center of plasma again
if par.determine_vector_potential
    TF_map.Z2=TF_map.Z2+Z_correction;   % Up the evaluation points by 0.07 since plasma is centered higher.
end
%% Loop over type of coils and nr. of coils
coil_name='TF_coil';
coils_struct=TF_coils;
for i=1:length(TF_coils)
    % Display which coil is calculated
    disp(['TF coil ',num2str(i)])
    
    % Coil geometry declaration
    BS.Nfilament=1;  % Number of filaments of a coil
    BS.filament.Gamma = [coils_struct(i).X,coils_struct(i).Y,coils_struct(i).Z]'; % Coordinates coil
    BS.filament.I = TF_map.I; % Current through coil
    BS.filament.dGamma = 1e-2; % Maximum size of coil element dl (otherwise 2D interpolation of coil coordinates)
    
    % Calculate B and store in TF coil map:
    [   TF_map.(coil_name)(i).BX,...
        TF_map.(coil_name)(i).BY,...
        TF_map.(coil_name)(i).BZ,...
        TF_map.(coil_name)(i).AX,...
        TF_map.(coil_name)(i).AY,...
        TF_map.(coil_name)(i).AZ] = ...
        BS_calc_B(BS,   TF_map.X,TF_map.Y,TF_map.Z,4);
    if par.determine_vector_potential
        [ ~,~,~,...
            AX2,...
            AY2,...
            AZ2] = BS_calc_B(BS,TF_map.X2,TF_map.Y2,TF_map.Z2,5); % Calculate only vector potential
        
        %% Make dA_dphi calculation
        AR2  =   TF_map.X2./TF_map.R2.*AX2...
            +TF_map.Y2./TF_map.R2.*AY2;
        Aphi2=  -TF_map.X2./TF_map.R2.*AY2...
            +TF_map.Y2./TF_map.R2.*AX2;
        
        TF_map.(coil_name)(i).dAR_dphi=...
            +1/12*AR2  (:,1)...
            -2/3 *AR2  (:,2)...
            +2/3 *AR2  (:,3)...
            -1/12*AR2  (:,4);
        TF_map.(coil_name)(i).dAZ_dphi=...
            +1/12*AZ2  (:,1)...
            -2/3 *AZ2  (:,2)...
            +2/3 *AZ2  (:,3)...
            -1/12*AZ2  (:,4);
        TF_map.(coil_name)(i).dAphi_dphi=...
            +1/12*Aphi2(:,1)...
            -2/3 *Aphi2(:,2)...
            +2/3 *Aphi2(:,3)...
            -1/12*Aphi2(:,4);
        
        %% Scale with size of angle
        TF_map.(coil_name)(i).dAR_dphi  =TF_map.(coil_name)(i).dAR_dphi  /par.delta_phi2;
        TF_map.(coil_name)(i).dAZ_dphi  =TF_map.(coil_name)(i).dAZ_dphi  /par.delta_phi2;
        TF_map.(coil_name)(i).dAphi_dphi=TF_map.(coil_name)(i).dAphi_dphi/par.delta_phi2;
    end
end
%% Remove obsolete fields
remove={'X2','Y2','Z2','R2'};
TF_map=remove_fields(TF_map,remove);
TF_map.Z=TF_map.Z-Z_correction;   % Reduce Z express the vertical distance from center of plasma again

disp('************************************************')
disp('DONE!')
disp('************************************************')
end