function [ ] = BS_AUG_RMP_main(PROCESS_NUMBER,N_PROCESS)
%%BS_AUG_RMP_main Determines the magnetic field of RMP field
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
        par.PROCESS_NUMBER=300;
        par.N_PROCESS=3000;
    case 1
        par.N_PROCESS=20;
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
    case 2
        par.N_PROCESS=str2double(N_PROCESS);
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
end

%% Parameters
RMP_map.I=4.8e3;        par.I=RMP_map.I;        % Coil current [kAt] (default),
par.determine_vector_potential=true;            % Determines derivatives in phi of vector potential. Needs auxillary points for nice finite difference calculation.
par.example=0;                                  % No plotting / example
par.coord_sys='toroidal';                           % Switch between coordinal systems (flux: theta,psi,phi OR toroidal: R,Z,phi)

%% Load Coils
[ RMP_coils] = BS_AUG_RMP_load_coils(1,1,0);

%% Define grid where to determine B
[RMP_map,Z_correction]=BS_AUG_get_grid(RMP_map);

%% Calculate the fields of each coil
[RMP_map] = BS_AUG_RMP_calc_coils_individual(par,RMP_map,RMP_coils,Z_correction);

remove={'index_link'};
par=remove_fields(par,remove);

%% Determine phi and R-components fields from original Cartesian
remove_obsolete={'BX','BY','AX','AY'}; % To remove obsolete fields
types=RMP_coils.coil_name;
for t=types
    coil_name=t{1};
    if ~isfield(RMP_map,coil_name)
%         warning(['Coils ',coil_name,' not present'])
        continue
    end
    for i=1:length(RMP_map.(coil_name))  % For each coil
        RMP_map.(coil_name)(i).Aphi  =-RMP_map.X./RMP_map.R.*RMP_map.(coil_name)(i).AY      +   RMP_map.Y./RMP_map.R.*RMP_map.(coil_name)(i).AX;
        RMP_map.(coil_name)(i).Bphi  =-RMP_map.X./RMP_map.R.*RMP_map.(coil_name)(i).BY      +   RMP_map.Y./RMP_map.R.*RMP_map.(coil_name)(i).BX;
        
        RMP_map.(coil_name)(i).AR    = RMP_map.X./RMP_map.R.*RMP_map.(coil_name)(i).AX      +   RMP_map.Y./RMP_map.R.*RMP_map.(coil_name)(i).AY;
        RMP_map.(coil_name)(i).BR    = RMP_map.X./RMP_map.R.*RMP_map.(coil_name)(i).BX      +   RMP_map.Y./RMP_map.R.*RMP_map.(coil_name)(i).BY;
        
        % Store neccesary fields in seperate variable and pre-allocate
        if i==1
            field.(coil_name)(length(RMP_map.(coil_name)))=remove_fields(RMP_map.(coil_name)(i),remove_obsolete); %#ok<*AGROW>
        end
        field.(coil_name)(i)=remove_fields(RMP_map.(coil_name)(i),remove_obsolete);
    end
end
%% Save output
if par.example==0
    save_name_individual=strcat('./output/BS_AUG_RMP_',par.coord_sys,'_individual_',datestr(now,'yyyy-mm-dd'),'_process',num2str(par.PROCESS_NUMBER),'.mat');
    save(save_name_individual,'-v7.3','field','par');
    disp(['Saved individual coil file: ',num2str(par.PROCESS_NUMBER),' of ',num2str(par.N_PROCESS)])
end


if par.example==0
    % OBSOLETE SAVING OF A MODE. THIS IS NOW DONE BASED ON INDIVIDUAL COIL
    % DATA (TO HAVE A VARIATY OF MODES AVAILABLE)
    
%     save_name_even2B=strcat('./output/BS_AUG_RMP_n=',num2str(tor_mode),'_even_',par.coord_sys,'Bfield_',datestr(now,'yyyy-mm-dd'),'_processmode',num2str(par.PROCESS_NUMBER),'.mat');
%     save_name_even2A=strcat('./output/BS_AUG_RMP_n=',num2str(tor_mode),'_even_',par.coord_sys,'Afield_',datestr(now,'yyyy-mm-dd'),'_processmode',num2str(par.PROCESS_NUMBER),'.mat');
%     fieldB=remove_fields(field,{'dAR_dphi','dAZ_dphi','dAphi_dphi','AR','AZ','Aphi'});
%     fieldA=remove_fields(field,{'BR','BZ','Bphi','AR','AZ','Aphi'});
%     save(save_name_even2B,'-v7.3','fieldB','par');
%     save(save_name_even2A,'-v7.3','fieldA','par');
%     disp(strcat('Saved n=',num2str(tor_mode),' file: ',num2str(par.PROCESS_NUMBER),' of ',num2str(par.N_PROCESS)))
    return
else
    %% EXAMPLES
    % Superposition of fields with tor mode and parity
    tor_mode=2;
    parity_name='even';
    RMP_coils=BS_AUG_RMP_add_signs(RMP_coils,tor_mode,parity_name);
    [ RMP_coils,field] = BS_AUG_RMP_apply_mode(RMP_map,RMP_coils);

    % Get into more easy format
    delete(findall(0,'type','figure','tag','BS_AUG_RMP_map_RMP_coils_figure'));
    drawnow;
    BX=NaN(par.size_total); BY=NaN(par.size_total); BZ=NaN(par.size_total);
    X_unper=NaN(par.size_total); Y_unper=NaN(par.size_total); Z_unper=NaN(par.size_total);
    BX(par.indexes)=field.BX(:);
    BY(par.indexes)=field.BY(:);
    BZ(par.indexes)=field.BZ(:);
    X_unper(par.indexes)=RMP_map.X(:); RMP_map.X=X_unper;
    Y_unper(par.indexes)=RMP_map.Y(:); RMP_map.Y=Y_unper;
    Z_unper(par.indexes)=RMP_map.Z(:); RMP_map.Z=Z_unper;
    for type=1:3
        for i=1:length(RMP_coils.(RMP_coils.coil_name{type}))
            BX_i=NaN(par.size_total); BY_i=NaN(par.size_total); BZ_i=NaN(par.size_total);
            BX_i(par.indexes)=RMP_map.(RMP_coils.coil_name{type})(i).BX(:);
            BY_i(par.indexes)=RMP_map.(RMP_coils.coil_name{type})(i).BY(:);
            BZ_i(par.indexes)=RMP_map.(RMP_coils.coil_name{type})(i).BZ(:);
            RMP_map.(RMP_coils.coil_name{type})(i).BX=BX_i;
            RMP_map.(RMP_coils.coil_name{type})(i).BY=BY_i;
            RMP_map.(RMP_coils.coil_name{type})(i).BZ=BZ_i;
        end
    end
    
%     BR=field.BR;
%     Bphi=field.Bphi;
    clear field
    
    % Calculate expressions which are easy
    B_tot=sqrt(BX.^2+BY.^2+BZ.^2);
    nabla_B=divergence(RMP_map.X,RMP_map.Y,RMP_map.Z,BX,BY,BZ); % Divergence
end
%% EXAMPLES
switch par.example
    case 0
        return
    %% B IN Y=0-PLANE
    case 1
        X=permute(X_unper,[3 2 1]); 	Y=permute(Y_unper,[3 2 1]);	Z=permute(Z_unper,[3 2 1]);
        Y_test=Y; 
        expr_y0=Y_test(2,2,:)==0;
        for type=1:3
            for i=1:length(RMP_coils.(RMP_coils.coil_name{type}))
                BX_i=RMP_map.(RMP_coils.coil_name{type})(i).BX;
                BY_i=RMP_map.(RMP_coils.coil_name{type})(i).BY;
                BZ_i=RMP_map.(RMP_coils.coil_name{type})(i).BZ;
                BX_i(isnan(BX_i))=0; BY_i(isnan(BY_i))=0; BZ_i(isnan(BZ_i))=0;
                nabla_B_i=divergence(X_unper,Y_unper,Z_unper,BX_i,BY_i,BZ_i); % Divergence
                BX_i=permute(BX_i,[3 2 1]); 	BY_i=permute(BY_i,[3 2 1]); 	BZ_i=permute(BZ_i,[3 2 1]);
                nabla_B_i=permute(nabla_B_i,[3 2 1]);
                [ha1,ha2,ha3,ha4]=make_figure([RMP_coils.coil_name{type},' coil ',num2str(i)]);
                
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BX_i(:,:,expr_y0),'parent',ha1,'displayname','$B_x$'), title(ha1,'$B_x$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BY_i(:,:,expr_y0),'parent',ha2,'displayname','$B_y$'), title(ha2,'$B_y$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),BZ_i(:,:,expr_y0),'parent',ha3,'displayname','$B_z$'), title(ha3,'$B_z$','interpreter','latex')
                contourf(X(:,:,expr_y0),Z(:,:,expr_y0),abs(nabla_B_i(:,:,expr_y0)),'parent',ha4,'displayname','$\nabla\vec{B}$'), title(ha4,'$|\nabla \vec{B}|$','interpreter','latex')
                colorbar('peer',ha1), colorbar('peer',ha2), colorbar('peer',ha3), colorbar('peer',ha4)
            end
        end
        %% Plot B of type of all
        for type=1:3
            if isempty(RMP_coils.(RMP_coils.coil_name{type})); continue; end
            BX_type=zeros(size(BX)); BY_type=zeros(size(BY)); BZ_type=zeros(size(BZ));
            
            for i=1:length(RMP_coils.(RMP_coils.coil_name{type}))
                BX_type=BX_type+RMP_coils.(RMP_coils.coil_name{type})(i).sign*RMP_map.(RMP_coils.coil_name{type})(i).BX;
                BY_type=BY_type+RMP_coils.(RMP_coils.coil_name{type})(i).sign*RMP_map.(RMP_coils.coil_name{type})(i).BY;
                BZ_type=BZ_type+RMP_coils.(RMP_coils.coil_name{type})(i).sign*RMP_map.(RMP_coils.coil_name{type})(i).BZ;
            end
            [ha1,ha2,ha3,ha4]=make_figure(['All ',RMP_coils.coil_name{type},' coils']);
            % Permute and plot
            nabla_B_type=divergence(RMP_map.X,RMP_map.Y,RMP_map.Z,BX_type,BY_type,BZ_type); % Divergence
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
        
        contourf(RMP_map.X(:,:,expr_z0),RMP_map.Y(:,:,expr_z0),abs(B_tot(:,:,expr_z0)),'parent',ha1,'displayname','$vec{B}$')
        
        colorbar('peer',ha1)
        xlabel(ha1,'$x$ [m]','interpreter','latex')
        ylabel(ha1,'$y$ [m]','interpreter','latex')
        title(ha1,'$|\vec{B}|$ in $z=0$-plane [T/m]','interpreter','latex')
        axis(ha1,'equal')
        
        plot_coils(ha1,RMP_coils)
        
        %% Plot X,Y-plane figure (top view torus) (divergence B)
        hf=figure('tag','BS_AUG_RMP_map_RMP_coils_figure');
        ha1=axes('parent',hf);
        hold(ha1,'on'), grid(ha1,'on')
        
        contourf(RMP_map.X(:,:,expr_z0),RMP_map.Y(:,:,expr_z0),abs(nabla_B(:,:,expr_z0)),'parent',ha1,'displayname','$|\nabla\vec{B}|$')
        
        colorbar('peer',ha1)
        xlabel(ha1,'$x$ [m]','interpreter','latex')
        ylabel(ha1,'$y$ [m]','interpreter','latex')
        title(ha1,'$|\nabla\vec{B}|$ in $z=0$-plane [T/m]','interpreter','latex')
        axis(ha1,'equal')
        plot_coils(ha1,RMP_coils)
    case 3
        X=permute(RMP_map.X,[3 2 1]); 	Y=permute(RMP_map.Y,[3 2 1]);	Z=permute(RMP_map.Z,[3 2 1]);
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
        
        plot_coils(ha1,RMP_coils)
        
        view(ha1,3)
        axis(ha1,'manual')
        %% Poloidal cross section (movie)
        
        for i=1:size(RMP_map.phi,3)
            % Slice psi-theta-phi coordinates
            RMP_map.X2=RMP_map.psi(:,:,i);
            RMP_map.Y2=RMP_map.theta(:,:,i);
            RMP_map.Z2=RMP_map.phi(:,:,i);
            % Plot this
            hslicer = slice(RMP_map.psi,RMP_map.theta,RMP_map.phi,B_tot,RMP_map.X2,RMP_map.Y2,RMP_map.Z2,'parent',ha1);
            set(hslicer,'EdgeColor','none')
            % Reset X-Y-Z-data
            set(hslicer,'XData',RMP_map.X(:,:,i),'YData',RMP_map.Y(:,:,i),'ZData',RMP_map.Z(:,:,i))
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
function plot_coils(ha1,RMP_coils)
Bu=RMP_coils.Bu;
Bl=RMP_coils.Bl;
A=RMP_coils.A;

for i=1:length(Bu)
    h=plot3(ha1,Bu(i).X,Bu(i).Y,Bu(i).Z,'r','displayname',['Bu coil ',num2str(i)]);
    if Bu(i).sign==1
        text(mean(Bu(i).X),mean(Bu(i).Y),mean(Bu(i).Z),'\textbf{+}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    elseif Bl(i).sign==-1
        text(mean(Bu(i).X),mean(Bu(i).Y),mean(Bu(i).Z),'\textbf{-}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    end
    set(h,'linewidth',3)
end
for i=1:length(Bl)
    h=plot3(ha1,Bl(i).X,Bl(i).Y,Bl(i).Z,'b','displayname',['Bl coil ',num2str(i)]);
    if Bl(i).sign==1
        text(mean(Bl(i).X),mean(Bl(i).Y),mean(Bl(i).Z),'\textbf{+}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    elseif Bl(i).sign==-1
        text(mean(Bl(i).X),mean(Bl(i).Y),mean(Bl(i).Z),'\textbf{-}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    end
    set(h,'linewidth',3)
end
for i=1:length(A)
    h=plot3(ha1,A(i).X,A(i).Y,A(i).Z,'k','displayname',['A coil ',num2str(i)]);
    if A(i).ref_coil==3
        set(h,'color','g')
    end
    if A(i).sign==1
        text(mean(A(i).X),mean(A(i).Y),mean(A(i).Z),'\textbf{+}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    elseif A(i).sign==-1
        text(mean(A(i).X),mean(A(i).Y),mean(A(i).Z),'\textbf{-}','color',get(h,'color'),'interpreter','latex','parent',ha1,'FontSize',20)
    end
    set(h,'linewidth',3)
end
end