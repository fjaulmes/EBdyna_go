function [RMP_map] = BS_AUG_RMP_calc_coils_individual(par,RMP_map,RMP_coils,Z_correction)
%BS_AUG_RMP_calc_coils_individual Summary of this function goes here
%   Detailed explanation goes here

%% Calculate B-field for each coil:
BS(1).Nfilament=1;

disp('************************************************')
disp('Starting standard B-field calculation of each coil:')
disp('************************************************')

%% Correction for offset magnetic axis
RMP_map.Z=RMP_map.Z+Z_correction;   % Reduce Z express the vertical distance from center of plasma again
if par.determine_vector_potential
    RMP_map.Z2=RMP_map.Z2+Z_correction;   % Up the evaluation points by 0.07 since plasme is centered higher.
end
%% Loop over type of coils and nr. of coils
for type=1:length(RMP_coils.coil_name)
    for i=1:length(RMP_coils.(RMP_coils.coil_name{type}))
        % Display which coil is calculated
        coil_name=RMP_coils.coil_name{type};
        coils_struct=RMP_coils.(RMP_coils.coil_name{type});
        disp([coil_name,' coil ',num2str(i)])
        
        % Coil geometry declaration
        BS.Nfilament=1;  % Number of filaments of a coil
        BS.filament.Gamma = [coils_struct(i).X,coils_struct(i).Y,coils_struct(i).Z]'; % Coordinates coil
        BS.filament.I = RMP_map.I; % Current through coil
        BS.filament.dGamma = 1e-3; % Maximum size of coil element dl (otherwise 2D interpolation of coil coordinates)
              
        % Calculate B and store in RMP_map:
        [   RMP_map.(coil_name)(i).BX,...
            RMP_map.(coil_name)(i).BY,...
            RMP_map.(coil_name)(i).BZ,...
            RMP_map.(coil_name)(i).AX,...
            RMP_map.(coil_name)(i).AY,...
            RMP_map.(coil_name)(i).AZ] = ...
                BS_calc_B(BS,   RMP_map.X,RMP_map.Y,RMP_map.Z,4);
        if par.determine_vector_potential
            [ ~,~,~,...
                AX2,...
                AY2,...
                AZ2] = BS_calc_B(BS,RMP_map.X2,RMP_map.Y2,RMP_map.Z2,5); % Calculate only vector potential
            
            %% Make dA_dphi calculation
            AR2  =   RMP_map.X2./RMP_map.R2.*AX2...
                    +RMP_map.Y2./RMP_map.R2.*AY2;
            Aphi2=  -RMP_map.X2./RMP_map.R2.*AY2...
                    +RMP_map.Y2./RMP_map.R2.*AX2;
            
            RMP_map.(coil_name)(i).dAR_dphi=...
                +1/12*AR2  (:,1)...
                -2/3 *AR2  (:,2)...
                +2/3 *AR2  (:,3)...
                -1/12*AR2  (:,4);
            RMP_map.(coil_name)(i).dAZ_dphi=...
                +1/12*AZ2  (:,1)...
                -2/3 *AZ2  (:,2)...
                +2/3 *AZ2  (:,3)...
                -1/12*AZ2  (:,4);
            RMP_map.(coil_name)(i).dAphi_dphi=...
                +1/12*Aphi2(:,1)...
                -2/3 *Aphi2(:,2)...
                +2/3 *Aphi2(:,3)...
                -1/12*Aphi2(:,4);
            
            %% Scale with size of angle
            RMP_map.(coil_name)(i).dAR_dphi  =RMP_map.(coil_name)(i).dAR_dphi  /par.delta_phi2;
            RMP_map.(coil_name)(i).dAZ_dphi  =RMP_map.(coil_name)(i).dAZ_dphi  /par.delta_phi2;
            RMP_map.(coil_name)(i).dAphi_dphi=RMP_map.(coil_name)(i).dAphi_dphi/par.delta_phi2;
        end
    end
end
%% Remove obsolete fields
remove={'X2','Y2','Z2','R2'};
RMP_map=remove_fields(RMP_map,remove);
RMP_map.Z=RMP_map.Z-Z_correction;   % Reduce Z express the vertical distance from center of plasma again

disp('************************************************')
disp('DONE!')
disp('************************************************')
end