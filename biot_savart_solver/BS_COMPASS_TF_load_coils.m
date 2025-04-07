function [ TF] = BS_COMPASS_TF_load_coils(INPUTFILE,NB_TF_COILS,REDUCE_DATA)
%BS_AUG_TFR_load_coils Load X,Y,Z TFR coils coordinates from data file
%   Return TF_coil data

%% Load files
if nargin==2
    REDUCE_DATA=1;
else
    if REDUCE_DATA>0
        warning(['reducing coil input data by a factor ' num2str(REDUCE_DATA)]);
    else
        REDUCE_DATA=1;
    end
end

if ~exist(INPUTFILE,'file')
    error('Map with COMPASS TF coil not found')
else
    % Pre-allocate
    TF(NB_TF_COILS).coil=[];
    % Load single coil data
    CUcoils=load(INPUTFILE);   
    
    
    if numel(fieldnames(CUcoils))==2
        CUcoils.X_coils=CUcoils.X_coils(1:REDUCE_DATA:end);
        CUcoils.Z_coils=CUcoils.Z_coils(1:REDUCE_DATA:end);
        % Please note the data consists of RZ coords for 1 TF coil
        
        coord_center=zeros(length(CUcoils.X_coils),3);
        coord_center(:,1)=CUcoils.X_coils;
        coord_center(:,3)=CUcoils.Z_coils;
        coord_center(:,2)=0;        % Define y=0 for basis coil, x -> R, z-> Z (so now coord_center has cylindrical coordinates)
        phi_rot_coils=linspace(0,2*pi,NB_TF_COILS+1); phi_rot_coils(end)=[];
        
        x_coils= coord_center(:,1)*cos(phi_rot_coils);
        y_coils=-coord_center(:,1)*sin(phi_rot_coils);  % Minus since of order / right-handed axis / phi runs in negative Y
        z_coils= coord_center(:,3)*ones(size(phi_rot_coils));
        
    
        % Put the data in coil struct
        for i=1:NB_TF_COILS
            TF(i).X=x_coils(:,i);
            TF(i).Y=y_coils(:,i);
            TF(i).Z=z_coils(:,i);
        end
    else
        CUcoils.X_coils=CUcoils.X_coils(1:REDUCE_DATA:end);
        CUcoils.Y_coils=CUcoils.Y_coils(1:REDUCE_DATA:end);
        CUcoils.Z_coils=CUcoils.Z_coils(1:REDUCE_DATA:end);
        disp('The data of the TF coils in provided as (X,Y,Z) and so a complete 3D wire is applied.');
        disp(['TF number is ' num2str(NB_TF_COILS)]);
        x_coils= CUcoils.X_coils';
        y_coils= CUcoils.Y_coils';
        z_coils= CUcoils.Z_coils';
        R_coils=hypot(x_coils,y_coils);
        phi_coils=acos(x_coils./R_coils);
        phi_rot_coils=linspace(0,2*pi,NB_TF_COILS+1); phi_rot_coils(end)=[];
        % Put the data in coil struct
        for i=1:NB_TF_COILS
            TF(i).X=R_coils.*cos(phi_coils+phi_rot_coils(i));
            TF(i).Y=R_coils.*sin(phi_coils+phi_rot_coils(i));
            TF(i).Z=z_coils;
        end
    end

end

%% UNIT TEST - plotting
if nargout==0
    delete(findall(0,'type','figure','tag',mfilename))
    hf=figure('tag',mfilename);
    ha=axes('parent',hf);
    hold(ha,'on'), grid(ha,'on'), view(ha,3)
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$y$ [m]','interpreter','latex')
    zlabel(ha,'$z$ [m]','interpreter','latex')
    title(ha,'UNIT TEST loading COMPASS TF coils','interpreter','latex')
    
    for i=1:NB_TF_COILS
        h=plot3(ha,TF(i).X,TF(i).Y,TF(i).Z,'k','displayname',['TF coil ',num2str(i)]);
        if i~=1
            hasbehavior(h,'legend',false)
        else
            set(h,'color','r')
        end
        
    end
    axis(ha,'equal')
    legend(ha,'show')
end

end