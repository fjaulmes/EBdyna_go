function [ ] = BS_AUG_TFR_main(PROCESS_NUMBER,N_PROCESS)
%%BS_AUG_TFR_main Determines the magnetic field of TFR field
% MAIN FILE
%   Script to determine field of coils in AUG. Can store field of each coil seperately in RMP_map-struct.
%   calls BS_AUG_RMP_load_coils to load in all coils from data in folder
%   define the number of psi, theta and phi points
%   Uses a standard current of 4.8kAt
global par
%#ok<*COLND>
%% PROCESS NUMBERS
switch nargin
    case 0
        par.PROCESS_NUMBER=1;
        par.N_PROCESS=3000;
    case 1
        par.N_PROCESS=20;
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
    case 2
        par.N_PROCESS=str2double(N_PROCESS);
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
end

%% Parameters
TF_map.I=0.3e6;   par.I=TF_map.I;       % Coil current [At] (default),
par.determine_vector_potential=true;    % Determines derivatives in phi of vector potential. Needs auxillary points for nice finite difference calculation.
par.example=0;                          % No plotting / example
par.coord_sys='toroidal';                   % Switch between coordinal systems
par.nr_coils=1;                         % Switch between 1 or 16 coils

%% Load Coils
[TF_coils] = BS_AUG_TF_load_coils;
TF_coils=TF_coils(1:par.nr_coils);

%% Define grid where to determine B
[TF_map,Z_correction]=BS_AUG_get_grid(TF_map);

%% Calculate the fields of each coil
[TF_map] = BS_AUG_TFR_calc_coils_individual(par,TF_map,TF_coils,Z_correction);

% Remove the link between A and B-positions / fields
remove_link={'index_link'};
par=remove_fields(par,remove_link);

%% Determine phi and R-components fields from original Cartesian
remove_obsolete={'BX','BY','AX','AY'}; % To remove obsolete fields
for i=1:length(TF_map.TF_coil)  % For each coil
    TF_map.TF_coil(i).Aphi  =-TF_map.X./TF_map.R.*TF_map.TF_coil(i).AY      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).AX;
    TF_map.TF_coil(i).Bphi  =-TF_map.X./TF_map.R.*TF_map.TF_coil(i).BY      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).BX;
    
    TF_map.TF_coil(i).AR    = TF_map.X./TF_map.R.*TF_map.TF_coil(i).AX      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).AY;
    TF_map.TF_coil(i).BR    = TF_map.X./TF_map.R.*TF_map.TF_coil(i).BX      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).BY;
    
    % Store neccesary fields in seperate variable and pre-allocate
    if i==1
        field(length(TF_map.TF_coil))=remove_fields(TF_map.TF_coil(i),remove_obsolete); %#ok<*AGROW>
    end
    field(i)=remove_fields(TF_map.TF_coil(i),remove_obsolete); 
end

%% Save output
if ~par.example
    save_name_individual=strcat('./output/BS_AUG_TF_',par.coord_sys,'_individual_',datestr(now,'yyyy-mm-dd'),'_process',num2str(par.PROCESS_NUMBER),'.mat');
    save(save_name_individual,'-v7.3','field','par');
    disp(['Saved individual coil file: ',num2str(par.PROCESS_NUMBER),' of ',num2str(par.N_PROCESS)])
    return
else
    %% EXAMPLES
    % Superposition of fields with tor mode and parity
    tor_mode=2;
    parity_name='even';
    TF_coils=BS_AUG_TFR_add_signs(TF_coils,tor_mode,parity_name);
    [ TF_coils,field] = BS_AUG_RMP_apply_mode(TF_map,TF_coils);

    % Get into more easy format
    delete(findall(0,'type','figure','tag','BS_AUG_RMP_map_RMP_coils_figure'));
    drawnow;
    BX=NaN(par.size_total); BY=NaN(par.size_total); BZ=NaN(par.size_total);
    X_unper=NaN(par.size_total); Y_unper=NaN(par.size_total); Z_unper=NaN(par.size_total);
    BX(par.indexes)=field.BX(:);
    BY(par.indexes)=field.BY(:);
    BZ(par.indexes)=field.BZ(:);
    X_unper(par.indexes)=TF_map.X(:); TF_map.X=X_unper;
    Y_unper(par.indexes)=TF_map.Y(:); TF_map.Y=Y_unper;
    Z_unper(par.indexes)=TF_map.Z(:); TF_map.Z=Z_unper;
    for type=1:3
        for i=1:length(TF_coils.(TF_coils.coil_name{type}))
            BX_i=NaN(par.size_total); BY_i=NaN(par.size_total); BZ_i=NaN(par.size_total);
            BX_i(par.indexes)=TF_map.(TF_coils.coil_name{type})(i).BX(:);
            BY_i(par.indexes)=TF_map.(TF_coils.coil_name{type})(i).BY(:);
            BZ_i(par.indexes)=TF_map.(TF_coils.coil_name{type})(i).BZ(:);
            TF_map.(TF_coils.coil_name{type})(i).BX=BX_i;
            TF_map.(TF_coils.coil_name{type})(i).BY=BY_i;
            TF_map.(TF_coils.coil_name{type})(i).BZ=BZ_i;
        end
    end
    
%     BR=field.BR;
%     Bphi=field.Bphi;
    clear field
    
    % Calculate expressions which are easy
    B_tot=sqrt(BX.^2+BY.^2+BZ.^2);
    nabla_B=divergence(TF_map.X,TF_map.Y,TF_map.Z,BX,BY,BZ); % Divergence
end
%% EXAMPLES (MADE for RMP - coils, not maintained)
switch par.example
    case 0
        return
    %% B IN Y=0-PLANE
    case 1
        X=permute(X_unper,[3 2 1]); 	Y=permute(Y_unper,[3 2 1]);	Z=permute(Z_unper,[3 2 1]);
        Y_test=Y; 
        expr_y0=Y_test(2,2,:)==0;
        for type=1:3
            for i=1:length(TF_coils.(TF_coils.coil_name{type}))
                BX_i=TF_map.(TF_coils.coil_name{type})(i).BX;
                BY_i=TF_map.(TF_coils.coil_name{type})(i).BY;
                BZ_i=TF_map.(TF_coils.coil_name{type})(i).BZ;
                BX_i(isnan(BX_i))=0; BY_i(isnan(BY_i))=0; BZ_i(isnan(BZ_i))=0;
                nabla_B_i=divergence(X_unper,Y_unper,Z_unper,BX_i,BY_i,BZ_i); % Divergence
                BX_i=permute(BX_i,[3 2 1]); 	BY_i=permute(BY_i,[3 2 1]); 	BZ_i=permute(BZ_i,[3 2 1]);
                nabla_B_i=permute(nabla_B_i,[3 2 1]);
                [ha1,ha2,ha3,ha4]=make_figure([TF_coils.coil_name{type},' coil ',num2str(i)]);
                
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BX_i(:,:,expr_y0),'parent',ha1,'displayname','$B_x$'), title(ha1,'$B_x$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BY_i(:,:,expr_y0),'parent',ha2,'displayname','$B_y$'), title(ha2,'$B_y$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BZ_i(:,:,expr_y0),'parent',ha3,'displayname','$B_z$'), title(ha3,'$B_z$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),abs(nabla_B_i(:,:,expr_y0)),'parent',ha4,'displayname','$\nabla\vec{B}$'), title(ha4,'$|\nabla \vec{B}|$','interpreter','latex')
                colorbar('peer',ha1), colorbar('peer',ha2), colorbar('peer',ha3), colorbar('peer',ha4)
            end
        end
        %% Plot B of type of all
        for type=1:3
            if isempty(TF_coils.(TF_coils.coil_name{type})); continue; end
            BX_type=zeros(size(BX)); BY_type=zeros(size(BY)); BZ_type=zeros(size(BZ));
            
            for i=1:length(TF_coils.(TF_coils.coil_name{type}))
                BX_type=BX_type+TF_coils.(TF_coils.coil_name{type})(i).sign*TF_map.(TF_coils.coil_name{type})(i).BX;
                BY_type=BY_type+TF_coils.(TF_coils.coil_name{type})(i).sign*TF_map.(TF_coils.coil_name{type})(i).BY;
                BZ_type=BZ_type+TF_coils.(TF_coils.coil_name{type})(i).sign*TF_map.(TF_coils.coil_name{type})(i).BZ;
            end
            [ha1,ha2,ha3,ha4]=make_figure(['All ',TF_coils.coil_name{type},' coils']);
            % Permute and plot
            nabla_B_type=divergence(TF_map.X,TF_map.Y,TF_map.Z,BX_type,BY_type,BZ_type); % Divergence
            nabla_B_type=permute(nabla_B_type,[3 2 1]);
            BX_type=permute(BX_type,[3 2 1]); BY_type=permute(BY_type,[3 2 1]); BZ_type=permute(BZ_type,[3 2 1]);
            contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BX_type(:,:,expr_y0),'parent',ha1,'displayname','$B_x$'), title(ha1,'$B_x$','interpreter','latex')
            contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BY_type(:,:,expr_y0),'parent',ha2,'displayname','$B_y$'), title(ha2,'$B_y$','interpreter','latex')
            contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BZ_type(:,:,expr_y0),'parent',ha3,'displayname','$B_z$'), title(ha3,'$B_z$','interpreter','latex')
            contourf(X(:,:,expr_y0),Z(:,:,expr_y0),abs(nabla_B_type(:,:,expr_y0)),'parent',ha4,'displayname','$\nabla\vec{B}$'), title(ha4,'$|\nabla\vec{B}|$','interpreter','latex')
            colorbar('peer',ha1), colorbar('peer',ha2), colorbar('peer',ha3), colorbar('peer',ha4)
        end
        %% Plot B of all coils together
        [ha1,ha2,ha3,ha4]=make_figure('All coils');
        % Permute and plot
        nabla_B=permute(nabla_B,[3 2 1]);
        BX_all=permute(BX,[3 2 1]); BY_all=permute(BY,[3 2 1]); BZ_all=permute(BZ,[3 2 1]);
        contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BX_all(:,:,expr_y0),'parent',ha1,'displayname','$B_x$'), title(ha1,'$B_x$','interpreter','latex')
        contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BY_all(:,:,expr_y0),'parent',ha2,'displayname','$B_y$'), title(ha2,'$B_y$','interpreter','latex')
        contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BZ_all(:,:,expr_y0),'parent',ha3,'displayname','$B_z$'), title(ha3,'$B_z$','interpreter','latex')
        contourf(X(:,:,expr_y0),Z(:,:,expr_y0),abs(nabla_B(:,:,expr_y0)),'parent',ha4,'displayname','$\nabla\vec{B}$'), title(ha4,'$|\nabla\vec{B}|$','interpreter','latex')
        colorbar('peer',ha1), colorbar('peer',ha2), colorbar('peer',ha3), colorbar('peer',ha4)
    case 2
        %% Plot X,Y-plane figure (top view torus) (B)
        hf=figure('tag','BS_AUG_RMP_map_RMP_coils_figure');
        ha1=axes('parent',hf);
        hold(ha1,'on'), grid(ha1,'on')
        
        contourf(TF_map.X(:,:,expr_z0),TF_map.Y(:,:,expr_z0),abs(B_tot(:,:,expr_z0)),'parent',ha1,'displayname','$vec{B}$')
        
        colorbar('peer',ha1)
        xlabel(ha1,'$x$ [m]','interpreter','latex')
        ylabel(ha1,'$y$ [m]','interpreter','latex')
        title(ha1,'$|\vec{B}|$ in $z=0$-plane [T/m]','interpreter','latex')
        axis(ha1,'equal')
        
        plot_coils(ha1,TF_coils)
        
        %% Plot X,Y-plane figure (top view torus) (divergence B)
        hf=figure('tag','BS_AUG_RMP_map_RMP_coils_figure');
        ha1=axes('parent',hf);
        hold(ha1,'on'), grid(ha1,'on')
        
        contourf(TF_map.X(:,:,expr_z0),TF_map.Y(:,:,expr_z0),abs(nabla_B(:,:,expr_z0)),'parent',ha1,'displayname','$|\nabla\vec{B}|$')
        
        colorbar('peer',ha1)
        xlabel(ha1,'$x$ [m]','interpreter','latex')
        ylabel(ha1,'$y$ [m]','interpreter','latex')
        title(ha1,'$|\nabla\vec{B}|$ in $z=0$-plane [T/m]','interpreter','latex')
        axis(ha1,'equal')
        plot_coils(ha1,TF_coils)
    case 3
        X=permute(TF_map.X,[3 2 1]); 	Y=permute(TF_map.Y,[3 2 1]);	Z=permute(TF_map.Z,[3 2 1]);
        %% Plot X,Y-plane figure (top view torus) (B)
        hf=figure('tag','BS_AUG_RMP_map_RMP_coils_figure');
        ha1=axes('parent',hf);
        hold(ha1,'on'), grid(ha1,'on')
        
        B_tot=permute(B_tot,[3 2 1]);
        contourf(X(:,:,expr_z0),Y(:,:,expr_z0),B_tot(:,:,expr_z0),'parent',ha1,'displayname','$\left|\vec{B}\right|$');
        c=get(gca,'children');
        axis(ha1,'equal')
        colorbar('peer',ha1)
        xlabel(ha1,'$x$ [m]','interpreter','latex')
        ylabel(ha1,'$y$ [m]','interpreter','latex')
        title(ha1,'$|\vec{B}|$ in $z=0$-plane [T/m]','interpreter','latex')
        %axis('equal')
        
        plot_coils(ha1,TF_coils)
        
        view(ha1,3)
        axis(ha1,'manual')
        %% Poloidal cross section (movie)
        
        for i=1:size(TF_map.phi,3)
            % Slice psi-theta-phi coordinates
            TF_map.X2=TF_map.psi(:,:,i);
            TF_map.Y2=TF_map.theta(:,:,i);
            TF_map.Z2=TF_map.phi(:,:,i);
            % Plot this
            hslicer = slice(TF_map.psi,TF_map.theta,TF_map.phi,B_tot,TF_map.X2,TF_map.Y2,TF_map.Z2,'parent',ha1);
            set(hslicer,'EdgeColor','none')
            % Reset X-Y-Z-data
            set(hslicer,'XData',TF_map.X(:,:,i),'YData',TF_map.Y(:,:,i),'ZData',TF_map.Z(:,:,i))
            col_b=colorbar('peer',ha1);
            
            drawnow;
            pause(0.1)
            if i~=2
                delete(hslicer)
            end
            delete(col_b)
        end
end
end

%% Function to create 2,2-subplot with certain name, X,Z-plane
function [ha1,ha2,ha3,ha4]=make_figure(name)
hf=figure('name',name,'tag','BS_AUG_RMP_map_RMP_coils_figure');
ha1=subplot(2,2,1,'parent',hf); axis(ha1,'equal')
ha2=subplot(2,2,2,'parent',hf); axis(ha2,'equal')
ha3=subplot(2,2,3,'parent',hf); axis(ha3,'equal')
ha4=subplot(2,2,4,'parent',hf); axis(ha4,'equal')
hold(ha1,'on'), hold(ha2,'on'),hold(ha3,'on'),hold(ha4,'on')
for ha=[ha1 ha2 ha3 ha4]
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$z$ [m]','interpreter','latex')
end
end

%% Function to add coil shape in figure
function plot_coils(ha1,TF)

for i=1:16
    h=plot3(ha1,TF(i).X,TF(i).Y,TF(i).Z,'k','displayname',['TF coil ',num2str(i)]);
    if i~=1
        hasbehavior(h,'legend',false)
    end
end
end