function [E,BR,BZ,Bphi] = B_interpolation(x,types)
global par dim maps time
%find_3D_Bfield Interpolates the 3D grid to find 3D field
%   Interpolates in 2D for flux coordinates, which are used to interpolate
%   for the 3D field(s).
%   Can use a variety of interpolation methods and indexes / slopes -
%   struct for 2D interpolation.
%   Support for different interpolation schemes:
% 0     -   self-written linear, with vector potential
% 1     -   self-written linear
% 2     -   interp linear               (ML standard)
% 3     -   ba_interp linear            (File Exchange)
% 4     -   griddedInterpolant linear   (makes use of lookup table)
% 5     -   interp cubic                (ML standard)
% 6     -   ba_interp cubic             (File Exchange)
% 7     -   griddedInterpolant cubic    (makes use of lookup table)
% 8     -   interp2 spline              (ML standard)
% 9     -   griddedInterpolant spline   (makes use of lookup table - RAM intensive!)
E=0;

%% UNIT TEST
if nargout == 0 && nargin ==0
    delete(findall(0,'type','figure','tag',mfilename));
    set(groot,'defaulttextinterpreter','latex')
    if par.APPLY_SAWTOOTH
        % UNIT TEST ST
        warning('UNIT TEST SAWTOOTH')
        delete(findall(0,'type','figure','tag',mfilename));
        
%         % Contour of LCFS
%         c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps.psi_norm_XZ',[dim.NB_PSI dim.NB_PSI]);
%         
%         % Make 2 figures
%         hf=figure('tag',mfilename,'name','Potential / flux','Windowstyle','normal');
%         hf2=figure('tag',mfilename,'name','Electric / magnetic fields','Windowstyle','normal');
%         
%         for i=1:6
%             ha1(i)=subplot(2,3,i,'parent',hf);
%             plot(ha1(i),c(1,2:end),(c(2,2:end)),'r','displayname','LCFS')
%             hold(ha1(i),'on');
%             axis(ha1(i),'equal')
%             
%             ha2(i)=subplot(2,3,i,'parent',hf2);
%             plot(ha2(i),c(1,2:end),(c(2,2:end)),'r','displayname','LCFS')
%             hold(ha2(i),'on');
%             axis(ha2(i),'equal')
%         end
%         linkaxes([ha1 ha2],'xy');
%         
%         % Set labels and titles
%         xlabel(ha1(4),'$R$','interpreter','latex'); xlabel(ha1(5),'$R$','interpreter','latex'); xlabel(ha1(6),'$R$','interpreter','latex');
%         ylabel(ha1(1),'$Z$','interpreter','latex'); ylabel(ha1(4),'$Z$','interpreter','latex');
%         xlabel(ha2(4),'$R$','interpreter','latex'); xlabel(ha2(5),'$R$','interpreter','latex'); xlabel(ha2(6),'$R$','interpreter','latex');
%         ylabel(ha2(1),'$Z$','interpreter','latex'); ylabel(ha2(4),'$Z$','interpreter','latex');
%         title(ha1(1),'$\psi (R,Z)$'); title(ha1(2),'$\psi(R+dR,Z)$'); title(ha1(3),'$\psi(R,Z+dZ)$')
%         title(ha1(4),'$\Phi (R,Z)$'); title(ha1(5),'$\Phi(R+dR,Z)$'); title(ha1(6),'$\Phi(R,Z+dZ)$')
%         title(ha2(1),'$B_R$'); title(ha2(2),'$B_Z$'); title(ha2(3),'$B_\varphi$')
%         title(ha2(4),'$E_R$'); title(ha2(5),'$E_Z$'); title(ha2(6),'$E_\varphi$')
%         
%         % Make 2D plane at phi=0;
%         X=linspace(1.5,1.9,1e2);
%         Z=linspace(-0.2,0.2,1e2);
%         [X2,Z2]=meshgrid(X,Z);
%         x=cat(2,X2(:),Z2(:),zeros(size(Z2(:))));
%         set(ha1(1),'Xlim',[min(X) max(X)],'Ylim',[min(Z) max(Z)])
%         
%         old_time=time;
%         time_movie=linspace(0,(par.st.t_reconnection+par.st.t_relaxation)*1.1,1e2);
%         for time=time_movie
%             if exist('h','var')
%                 delete(h);
%             end
%             %             error('Add Phi and psi_star to output arguments AND set expr_st in B_ST to true(size(x,1));')
%             [E,B,~,~]=B_interpolation(x,'ST_2D');  %PLEASE SET EXPR_ST TO TRUE FOR EACH PARTICLE TO PRODUCE PROPER IMAGE
%             if isempty(Phi)
%                 Phi=zeros(size(psi_star));
%             end
%             if numel(E)==1
%                 E=zeros(size(B));
%             end
%             
%             % Figure 1
%             for pos=1:3
%                 Phi_plot=reshape(Phi(:,:,:,pos),size(X2));
%                 psi_star_plot=reshape(psi_star(:,:,:,pos),size(X2));
%                 [h(1,pos)]=imagesc(X,Z,psi_star_plot,'parent',ha1(pos));
%                 [h(1,pos+3)]=imagesc(X,Z,Phi_plot,'parent',ha1(pos+3));
%             end
%             % Figure 2
%             BR  =reshape(B(:,1),size(X2));
%             BZ  =reshape(B(:,2),size(X2));
%             Bphi=reshape(B(:,3),size(X2));
%             ER  =reshape(E(:,1),size(X2));
%             EZ  =reshape(E(:,2),size(X2));
%             Ephi=reshape(E(:,3),size(X2));
%             
%             [h(2,1)]=imagesc(X,Z,BR  ,'parent',ha2(1));
%             [h(2,2)]=imagesc(X,Z,BZ  ,'parent',ha2(2));
%             [h(2,3)]=imagesc(X,Z,Bphi,'parent',ha2(3));
%             [h(2,4)]=imagesc(X,Z,ER  ,'parent',ha2(4));
%             [h(2,5)]=imagesc(X,Z,EZ  ,'parent',ha2(5));
%             [h(2,6)]=imagesc(X,Z,Ephi,'parent',ha2(6));
%             
%             for i=1:6
%                 colorbar('peer',ha1(i));
%                 colormap(ha1(i),'jet');
%                 colorbar('peer',ha2(i));
%                 colormap(ha2(i),'jet');
%             end
%             h(2,7) =quiver(X2(10:10:end,10:10:end),Z2(10:10:end,10:10:end),BR(10:10:end,10:10:end),BZ(10:10:end,10:10:end),'parent',ha2(1),'color','k');
%             h(2,8) =quiver(X2(10:10:end,10:10:end),Z2(10:10:end,10:10:end),BR(10:10:end,10:10:end),BZ(10:10:end,10:10:end),'parent',ha2(2),'color','k');
%             h(2,9) =quiver(X2(10:10:end,10:10:end),Z2(10:10:end,10:10:end),ER(10:10:end,10:10:end),EZ(10:10:end,10:10:end),'parent',ha1(4),'color','k');
%             h(2,10)=quiver(X2(10:10:end,10:10:end),Z2(10:10:end,10:10:end),ER(10:10:end,10:10:end),EZ(10:10:end,10:10:end),'parent',ha1(5),'color','k');
%             h(2,11)=quiver(X2(10:10:end,10:10:end),Z2(10:10:end,10:10:end),ER(10:10:end,10:10:end),EZ(10:10:end,10:10:end),'parent',ha1(6),'color','k');
%             
%             drawnow;
%         end
%         time=old_time;
%         warning('UNIT TEST COMPLETED')
    else
%         hf=figure('tag',mfilename);
%         
%         ha1(1,1)=subplot(3,2,1,'parent',hf); ha1(1,2)=subplot(3,2,2,'parent',hf);
%         ha1(2,1)=subplot(3,2,3,'parent',hf); ha1(2,2)=subplot(3,2,4,'parent',hf);
%         ha1(3,1)=subplot(3,2,5,'parent',hf); ha1(3,2)=subplot(3,2,6,'parent',hf);
%         
%         % Make 2D plane at phi=0;
%         X=dim.scale_X(5:end-5)+dim.R0+dim.DX/2;
%         Z=dim.scale_Z(5:end-5)+dim.DZ/2;
%         [X,Z]=ndgrid(X,Z);
%         x=cat(2,X(:),Z(:),zeros(size(Z(:))));
%         
%         [~,B_2D_BR,B_2D_BZ,B_2D_Bphi]=B_interpolation(x,'2D');
%         [~,B]=B_interpolation(x,'2D');
%         if ~isequaln(B,cat(2,B_2D_BR,B_2D_BZ,B_2D_Bphi))
%             error('Problem in conversions from size B-field')
%         end
%         B_2D_BR     =reshape(B_2D_BR    ,[length(dim.scale_X(5:end-5)),length(dim.scale_Z(5:end-5))]);
%         B_2D_BZ     =reshape(B_2D_BZ    ,[length(dim.scale_X(5:end-5)),length(dim.scale_Z(5:end-5))]);
%         B_2D_Bphi   =reshape(B_2D_Bphi  ,[length(dim.scale_X(5:end-5)),length(dim.scale_Z(5:end-5))]);
%         
%         imagesc(B_2D_BR'           ,'parent',ha1(1,2));
%         imagesc(B_2D_BZ'           ,'parent',ha1(2,2));
%         imagesc(B_2D_Bphi'         ,'parent',ha1(3,2));
%         
%         % Control
%         if any(par.interp_scheme==[3 6])
%             imagesc(maps(1).B_2D(:,:,1)','parent',ha1(1,1));
%             imagesc(maps(1).B_2D(:,:,2)','parent',ha1(2,1));
%             imagesc(maps(1).B_2D(:,:,3)','parent',ha1(3,1));
%         else
%             imagesc(maps(1).B_2D.BR'    ,'parent',ha1(1,1));
%             imagesc(maps(1).B_2D.BZ'    ,'parent',ha1(2,1));
%             imagesc(maps(1).B_2D.Bphi'  ,'parent',ha1(3,1));
%         end
%         
%         linkprop(ha1(1,:),'Clim');
%         linkprop(ha1(2,:),'Clim');
%         linkprop(ha1(3,:),'Clim');
%         linkprop(ha1(1,:),'Xdir'); linkprop(ha1(1,:),'Ydir');
%         linkprop(ha1(2,:),'Xdir'); linkprop(ha1(2,:),'Ydir');
%         linkprop(ha1(3,:),'Xdir'); linkprop(ha1(3,:),'Ydir');
%         colorbar('peer',ha1(1,1)); colorbar('peer',ha1(1,2));
%         colorbar('peer',ha1(2,1)); colorbar('peer',ha1(2,2));
%         colorbar('peer',ha1(3,1)); colorbar('peer',ha1(3,2));
%         
%         % Make similar figures for 3D interpolation at Z=0
%         if par.APPLY_3D
%             figure('tag',mfilename);
%             ha1(1,1)=subplot(3,2,1); ha1(1,2)=subplot(3,2,2);
%             ha1(2,1)=subplot(3,2,3); ha1(2,2)=subplot(3,2,4);
%             ha1(3,1)=subplot(3,2,5); ha1(3,2)=subplot(3,2,6);
%             if strcmp(par.coord_syst,'flux')
%                 m=load([par.paths.DATA_FOLDER,'flux_geometry.mat'],'X_PR_map','Z_PR_map');
%                 X=m.X_PR_map(1:513,1:513)+dim.R0;
%                 Z=m.Z_PR_map(1:513,1:513);
%             else
%                 X=dim.scale_X;
%                 Z=dim.scale_Z;
%                 [X,Z]=ndgrid(X,Z);
%             end
%             x=cat(2,X(:),Z(:),zeros(size(Z(:))));
%             
%             [~,B_3D_BR,B_3D_BZ,B_3D_Bphi]=B_interpolation(x,'3D');
%             B_3D_BR     =reshape(B_3D_BR    ,[513 513]);
%             B_3D_BZ     =reshape(B_3D_BZ    ,[513 513]);
%             B_3D_Bphi   =reshape(B_3D_Bphi  ,[513 513]);
%             
%             imagesc(B_3D_BR'    ,'parent',ha1(1,2));
%             imagesc(B_3D_BZ'  	,'parent',ha1(2,2));
%             imagesc(B_3D_Bphi'  ,'parent',ha1(3,2));
%             
%             if any(par.interp_scheme==[3 6])
%                 imagesc(maps(1).B_3D(:,:,1,1)'	,'parent',ha1(1,1));
%                 imagesc(maps(1).B_3D(:,:,1,2)'    ,'parent',ha1(2,1));
%                 imagesc(maps(1).B_3D(:,:,1,3)'  ,'parent',ha1(3,1));
%             else
%                 imagesc(maps(1).n3D.BR(:,:,1)'	,'parent',ha1(1,1));
%                 imagesc(maps(1).n3D.BZ(:,:,1)'    ,'parent',ha1(2,1));
%                 imagesc(maps(1).n3D.Bphi(:,:,1)'  ,'parent',ha1(3,1));
%             end
%             
%             linkprop(ha1(1,:),'Clim');
%             linkprop(ha1(2,:),'Clim');
%             linkprop(ha1(3,:),'Clim');
%             linkprop(ha1(1,:),'Xdir'); linkprop(ha1(1,:),'Ydir');
%             linkprop(ha1(2,:),'Xdir'); linkprop(ha1(2,:),'Ydir');
%             linkprop(ha1(3,:),'Xdir'); linkprop(ha1(3,:),'Ydir');
%             colorbar('peer',ha1(1,1)); colorbar('peer',ha1(1,2));
%             colorbar('peer',ha1(2,1)); colorbar('peer',ha1(2,2));
%             colorbar('peer',ha1(3,1)); colorbar('peer',ha1(3,2));
%         end
    end
    return
end

%% Switch the 'types', defining the intepolation field
if nargin<2
    if par.APPLY_SAWTOOTH && ~par.APPLY_3D
        % Sawtooth, but only in 2D
        types='ST_2D';
    elseif par.APPLY_SAWTOOTH
        % Sawtooth, and an extra interpolation for the 3D-field VACUUM without equilibrium!
        types='ST_3D';
    elseif par.APPLY_3D
        % The 3D simulation WITH equilibrium (pre-superpositioned in toroidal coordinates)
        types=par.coord_syst;
        % removed E=0 to consider centrifugal effects (Fc_field)
        E=x*0;
    else
        % The initial 2D equilibrium
        types='2D';
        % removed E=0 to consider centrifugal effects (Fc_field)
        E=x*0;
    end
end

%% Execute the 2D and/or 3D interpolation functions
switch types
    case '2D'
        [BR,BZ,Bphi]=B_2D(x);
    case 'flux'
        [BR,BZ,Bphi,X_ind,Z_ind]    =B_2D(x);
        [BR2,BZ2,Bphi2]             =B_3D_flux(x,X_ind,Z_ind);
        %Superimpose 2D and 3D fields
        BR  =   BR  +   BR2;
        BZ  =   BZ  +   BZ2;
        Bphi=   Bphi+   Bphi2;
    case 'toroidal'
        if par.superimpose_2D_3D
            [BR,BZ,Bphi]    =B_3D_toroidal(x);
        else
            [BR,BZ,Bphi,X_ind,Z_ind] = B_2D(x);
            [BR_3D,BZ_3D,Bphi_3D]    = B_3D_toroidal(x,X_ind,Z_ind);
            if isempty(BZ)
                BR=BR+BR_3D;
            else
                BR=BR+BR_3D;
                BZ=BZ+BZ_3D;
                Bphi=Bphi+Bphi_3D;
            end
        end
    case 'ST_2D'
        [BR,E]=B_ST(x);
        BZ=[]; Bphi=[];
    case 'ST_3D'
        [BR,E]    =B_ST(x);
        [B3D]     =B_3D_toroidal(x);
        
        % Combining with multiple time points (evaluation) requires a permutation:
        if size(x,3)==1
            B3D=permute(B3D,[1 3 2]);
        else
            B3D=permute(B3D,[1 4 3 2]);
        end
        
        BR=BR+B3D;
        BZ=[]; Bphi=[];
end

%% Adjust size
% Make the ordering: || particle | direction | time stamp || (if not already done by ST)
if any(par.interp_scheme==[3 6]) && ~any(strcmp(types,{'ST_2D','ST_3D'}))
    if size(x,3)==1
        BR  = permute(BR,[1 3 2]);
    else
        BR  =  permute(BR  ,[1 4 3 2]);
    end
end

% Convert from or to matrix with all components, or 3 vectors with each a
% component
if any(par.interp_scheme==[3 6]) && (nargout==4 || strcmp(par.scheme,'FabienB'))
    BZ  =BR(:,2,:);
    Bphi=BR(:,3,:);
    BR  =BR(:,1,:);
elseif ~any(par.interp_scheme==[3 6]) && (strcmp(par.scheme,'BORIS') || nargout<4)
    BR=cat(2,BR,BZ,Bphi);
end

end


%% B (2D)
function [BR,BZ,Bphi,X_ind,Z_ind]        =B_2D(x)
global par maps dim
persistent size_2D  F G H LAST_SCHEME

% Indexes 2D
X_ind=((x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;

switch par.interp_scheme
    case 0 	%% SELF-WRITTEN LINEAR DIVERGENCE FREE / POTENTIAL BASED
        R_inv=1./x(:,1,:);
        [BR,BZ,Bphi] = get_B_div_free_2D(X_ind,Z_ind,R_inv);
    case 1  %% SELF-WRITTEN LINEAR
        if isempty(size_2D)
            size_2D=[dim.size_X,dim.size_Z];
        end
        [X_ind,Z_ind] = interp_index_list_2D (size_2D,X_ind,Z_ind);     % Return the indexes / slopes as reps. X_ind and Z_ind
        if size(x,3)==1
            BR  =interp2_XZ(X_ind,Z_ind,maps(1).B_2D.BR);
            BZ  =interp2_XZ(X_ind,Z_ind,maps(1).B_2D.BZ);
            Bphi=interp2_XZ(X_ind,Z_ind,maps(1).B_2D.Bphi);
        else
            BR=NaN(size(Z_ind.x)); BZ=NaN(size(Z_ind.x)); Bphi=NaN(size(Z_ind.x));
            expr_c=isfinite(Z_ind.x);
            BR(expr_c)  =interp2_XZ(X_ind,Z_ind,maps(1).B_2D.BR  ,expr_c(:));
            BZ(expr_c)  =interp2_XZ(X_ind,Z_ind,maps(1).B_2D.BZ  ,expr_c(:));
            Bphi(expr_c)=interp2_XZ(X_ind,Z_ind,maps(1).B_2D.Bphi,expr_c(:));
        end
    case 2  %% INTERP2 LINEAR
        BR  =interp2(maps(1).B_2D.BR     ,Z_ind,X_ind,'*linear');
        BZ  =interp2(maps(1).B_2D.BZ     ,Z_ind,X_ind,'*linear');
        Bphi=interp2(maps(1).B_2D.Bphi   ,Z_ind,X_ind,'*linear');
    case 3  %% BA_INTERP2 LINEAR
        BR=ba_interp2(maps(1).B_2D       ,Z_ind,X_ind,'linear');
        BZ  =[];
        Bphi=[];
    case 4  %GRIDDEDINTERPOLANT LINEAR
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).B_2D.BR    ,'linear');
            G=griddedInterpolant(maps(1).B_2D.BZ    ,'linear');
            H=griddedInterpolant(maps(1).B_2D.Bphi  ,'linear');
        end
        BR  =F(X_ind,Z_ind);
        BZ  =G(X_ind,Z_ind);
        Bphi=H(X_ind,Z_ind);
    case 5  %% INTERP2 CUBIC
        BR  =interp2(maps(1).B_2D.BR     ,Z_ind,X_ind,'*cubic');
        BZ  =interp2(maps(1).B_2D.BZ     ,Z_ind,X_ind,'*cubic');
        Bphi=interp2(maps(1).B_2D.Bphi   ,Z_ind,X_ind,'*cubic');
    case 6  %% BA_INTERP2 CUBIC
        BR=ba_interp2(maps(1).B_2D       ,Z_ind,X_ind,'cubic');
        BZ  =[];
        Bphi=[];
    case 7  %GRIDDEDINTERPOLANT CUBIC
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).B_2D.BR    ,'cubic');
            G=griddedInterpolant(maps(1).B_2D.BZ    ,'cubic');
            H=griddedInterpolant(maps(1).B_2D.Bphi  ,'cubic');
        end
        BR  =F(X_ind,Z_ind);
        BZ  =G(X_ind,Z_ind);
        Bphi=H(X_ind,Z_ind);
    case 8  %% INTERP2 SPLINE
        BR  =interp2(maps(1).B_2D.BR     ,Z_ind,X_ind,'*spline');
        BZ  =interp2(maps(1).B_2D.BZ     ,Z_ind,X_ind,'*spline');
        Bphi=interp2(maps(1).B_2D.Bphi   ,Z_ind,X_ind,'*spline');
    case 9  %GRIDDEDINTERPOLANT SPLINE
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).B_2D.BR    ,'spline');
            G=griddedInterpolant(maps(1).B_2D.BZ    ,'spline');
            H=griddedInterpolant(maps(1).B_2D.Bphi  ,'spline');
        end
        BR  =F(X_ind,Z_ind);
        BZ  =G(X_ind,Z_ind);
        Bphi=H(X_ind,Z_ind);
    otherwise
        error('Interpolation scheme not written yet')
end
LAST_SCHEME=par.interp_scheme;

    function [BR,BZ,Bphi] = get_B_div_free_2D(X_ind,Z_ind,R_inv)
        % The index in the matrices / grid of the point 000 (i.e. where all
        % Beta / Gamma / Alpha's are stored
        X_ind_floor=floor(X_ind);
        Z_ind_floor=floor(Z_ind);
        
        % The local coordinates (from 0 to DX)
        dr  =(X_ind  -X_ind_floor)  ;
        dz  =(Z_ind  -Z_ind_floor)  ;
        
        % Fast sub2ind
        index = X_ind_floor+(Z_ind_floor-1)*(dim.size_X-1);
        index(isnan(index))=1;  % For ejected particles, retrieve the first element in matrices, which should be NAN!)
        
        % Calculation steps
        BR   = (maps(1).n2D.Beta_r  (index) + dr .* maps(1).n2D.Alpha_r  (index) ).*R_inv;
        BZ   = (maps(1).n2D.Beta_z  (index) + dz .* maps(1).n2D.Alpha_z  (index) ).*R_inv;
        Bphi =  maps(1).n2D.Beta_phi(index) + dr .* maps(1).n2D.Beta_phi_r(index) + dz .* (maps(1).n2D.Beta_phi_z(index) + dr .* maps(1).n2D.Beta_phi_r_z(index));
    end
end


%% B (SAWTOOTH)
function [B,E]        =B_ST(x)
global par maps dim time
persistent F G H CHI LAST_SCHEME

% The end evaluation is done with a for-loop, since the time is of
% importance in the ST simulation!
if size(x,3)>1
    B=zeros(size(x));
    E=zeros(size(x));
    old_time=time;
    time=par.dt; % First one to produce a figure (if configured)
    [~,~]=B_ST(x(:,:,1));
    for i=1:size(x,3)
        time=par.time_scale(i);
        [B(:,:,i),E(:,:,i)]=B_ST(x(:,:,i));
    end
    time=old_time;
    return
end

% Devide particles in 2 (sawtooth and equilibrium)
% 2D indexes for psi
X_ind= (x(:,1)-dim.R0)*dim.DX_inv +dim.mid_Xzero;
Z_ind=  x(:,2)        *dim.DZ_inv +dim.mid_Z;

psi_ind = ba_interp2(maps(1).psi_norm_XZ,Z_ind,X_ind,'linear');
% expr_st= true(size(psi_ind));
expr_st= psi_ind<1.1*dim.st.ind_rmix;    % Logical for particles within sawtooth region

% 2d indexes to determine theta
X_ind_dR= (x(expr_st,1)-dim.R0+par.st.DR)*dim.DX_inv +dim.mid_Xzero;
Z_ind_dZ= (x(expr_st,2)       +par.st.DZ)*dim.DZ_inv +dim.mid_Z;

X_ind_vec=cat(4,X_ind(expr_st,:),X_ind_dR        ,X_ind(expr_st,:));
Z_ind_vec=cat(4,Z_ind(expr_st,:),Z_ind(expr_st,:),Z_ind_dZ);

theta_vec= interpolate_theta_XZ(X_ind_vec,Z_ind_vec);

%% Equilibrium field
switch par.interp_scheme
    case 3  %% BA_INTERP2 LINEAR
        B=permute(ba_interp2(maps(1).B_2D       ,Z_ind,X_ind,'linear'),[1 3 2]); % NOTE: PERMUTATION SINCE SQUEEZE CANNOT DEAL WITH 1 PARTICLE AND E NEEDS TO HAVE PROPER SIZE. ONLY VALID FOR 1 TIMESTEP AT A TIME
        
        % Determine toroidal flux for particles within sawtooth region
        chi_vec= ba_interp2(maps(1).chi_XZ,Z_ind_vec,X_ind_vec,'linear');
    case 4  %GRIDDEDINTERPOLANT LINEAR
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).B_2D.BR    ,'linear');
            G=griddedInterpolant(maps(1).B_2D.BZ    ,'linear');
            H=griddedInterpolant(maps(1).B_2D.Bphi  ,'linear');
            CHI = griddedInterpolant(maps(1).chi_XZ ,'linear');
        end
        B   =F(X_ind,Z_ind);
        BZ  =G(X_ind,Z_ind);
        Bphi=H(X_ind,Z_ind);
        B = cat(2,B,BZ,Bphi);
        
        chi    = CHI(X_ind(expr_st,:) ,Z_ind(expr_st,:));
        chi_dR = CHI(X_ind_dR           ,Z_ind(expr_st,:));
        chi_dZ = CHI(X_ind(expr_st,:) ,Z_ind_dZ);
        
        chi_vec=cat(4,chi,chi_dR,chi_dZ);
    otherwise
        error('Interpolation scheme supported in sawtooth simulation')
end
E=zeros(size(B));

% Early return if sawtooth hasn't started yet (output parameters for if script has been adjusted for UNIT TEST)
if time==0; Phi_vec=zeros(size(chi_vec));psi_vec=zeros(size(chi_vec)); return; end

%% Sawtooth potential and flux -functions
% Determine r-coordinate for sawtooth particles
r_vec= sqrt(2*abs(chi_vec));

% Determine x and y-coordinates of particle in question
x_vec = r_vec .* cos(bsxfun(@minus,theta_vec,x(expr_st,3)));
y_vec = r_vec .* sin(bsxfun(@minus,theta_vec,x(expr_st,3)));

% Sawtooth parameters (time dependent)
st = st_parameters;
% Potential and psi^* (time + particle / space dependent)
[Phi_vec,psi_star_vec] = get_potential_flux(st,x_vec,y_vec,r_vec.^2);

% Convert / determine parameters for E and B
% helical angle calculation for debugging
% theta_st=squeeze(theta_vec(:,:,:,1));
% phi_st=mod(squeeze(x(expr_st,3)),2*pi);
% omega=theta_st-phi_st;
% PHI_POP=boolean((abs(phi_st)<0.3));
% SMALL_OMEGA=boolean((abs(omega)<0.3).*(abs(theta_st)<1.5));
% x_st=squeeze(x(expr_st,1));
% z_st=squeeze(x(expr_st,2));

psi_vec=psi_star_vec+chi_vec;
R_inv=1./x(expr_st,1,:);

%% Sawtooth fields
% Poloidal B-field from flux. Do not touch toroidal field.
if time<dim.st.time(par.st.stable_time)
    % Make a linear interpolation up to the first time stap
    B2=B; % Temporarely copy the equilibrium field
%     Bstar=B;
%     BH=B;
    B2(expr_st,1)=(-diff(psi_vec(:,:,:,[1 3]),1,4)).*R_inv/par.st.DZ;
    B2(expr_st,2)=( diff(psi_vec(:,:,:,1:2)  ,1,4)).*R_inv/par.st.DR;
%     BH(expr_st,1)=(-diff(chi_vec(:,:,:,[1 3]),1,4)).*R_inv/par.st.DZ;
%     BH(expr_st,2)=( diff(chi_vec(:,:,:,1:2)  ,1,4)).*R_inv/par.st.DR;
%     Bstar(expr_st,1)=(-diff(psi_star_vec(:,:,:,[1 3]),1,4)).*R_inv/par.st.DZ;
%     Bstar(expr_st,2)=( diff(psi_star_vec(:,:,:,1:2)  ,1,4)).*R_inv/par.st.DR;
%     B2(expr_st,1)=BH(expr_st,1)+Bstar(expr_st,1);
%     B2(expr_st,2)=BH(expr_st,2)+Bstar(expr_st,2);
%     BstarX=squeeze(Bstar(expr_st,1));
%     BstarZ=squeeze(Bstar(expr_st,2));
%     BX_st=squeeze(B(expr_st,1));
%     BZ_st=squeeze(B(expr_st,2));
%     BHX_st=squeeze(BH(expr_st,1));
%     BHZ_st=squeeze(BH(expr_st,2));
%     B2X_st=squeeze(B2(expr_st,1));
%     B2Z_st=squeeze(B2(expr_st,2));    
    
    % Poloidal E-field from potential. Toroidal field from dot product.
    E(expr_st,1)=dim.st.sign_psi_pol*(-diff(Phi_vec(:,:,:,1:2)  ,1,4)/par.st.DR);
    E(expr_st,2)=dim.st.sign_psi_pol*(-diff(Phi_vec(:,:,:,[1 3]),1,4)/par.st.DZ);
    E(expr_st,3)=-dot(E(expr_st,1:2),B2(expr_st,1:2),2)./B2(expr_st,3);
%     EX_st=squeeze(E(expr_st,1));
%     EZ_st=squeeze(E(expr_st,2));  
%     Ephi_st=squeeze(E(expr_st,3)); 
    
    t_slope=time/dim.st.time(par.st.stable_time);
    B(expr_st,:)=(1-t_slope)*B(expr_st,:)+t_slope*B2(expr_st,:);
    E(expr_st,:)=t_slope*E(expr_st,:);
elseif time>dim.st.time(end)
    % Use the sawtooth function from the first timestamp onwards.
    B(expr_st,1)=(-diff(psi_vec(:,:,:,[1 3]),1,4)).*R_inv/par.st.DZ;
    B(expr_st,2)=( diff(psi_vec(:,:,:,1:2)  ,1,4)).*R_inv/par.st.DR;
else
    % Use the sawtooth function from the first timestamp onwards.
    B(expr_st,1)=(-diff(psi_vec(:,:,:,[1 3]),1,4)).*R_inv/par.st.DZ;
    B(expr_st,2)=( diff(psi_vec(:,:,:,1:2)  ,1,4)).*R_inv/par.st.DR;
    
    % Poloidal E-field from potential. Toroidal field from dot product.
    E(expr_st,1)=dim.st.sign_psi_pol*(-diff(Phi_vec(:,:,:,1:2)  ,1,4)/par.st.DR);
    E(expr_st,2)=dim.st.sign_psi_pol*(-diff(Phi_vec(:,:,:,[1 3]),1,4)/par.st.DZ);
    E(expr_st,3)=-dot(E(expr_st,1:2),B(expr_st,1:2),2)./B(expr_st,3);
%     EX_st=squeeze(E(expr_st,1));
%     EZ_st=squeeze(E(expr_st,2));  
%     Ephi_st=squeeze(E(expr_st,3)); 
%     [xq,yq] = meshgrid(2.4:0.01:3.3, -0.5:0.01:0.5);
%     vq=griddata(x_st(PHI_POP),z_st(PHI_POP),squeeze(Phi_vec(PHI_POP,:,:,1)),xq,yq);
end
LAST_SCHEME=par.interp_scheme;
%% Sawtooth function 1 (time dependent parameters)
    function [st] = st_parameters
        persistent R2 R2_dot K_M_KR_dot K
        persistent ha h plotted
        if isempty(R2) || time==par.dt
            R2          =griddedInterpolant(dim.st.time,dim.st.r2        ,'linear');
            R2_dot      =griddedInterpolant(dim.st.time,dim.st.r2_dot    ,'linear');
            K           =griddedInterpolant(dim.st.time,dim.st.k         ,'linear');
            K_M_KR_dot  =griddedInterpolant(dim.st.time,dim.st.k_m_kr_dot,'linear');
            if par.st.show_parameters
                delete(findall(0,'type','figure','tag',mfilename));
                hf=figure('tag',mfilename,'windowstyle','normal','units','normalized','position',[0 0 1 1]);
                ha=NaN(9,1);
                plotted={'r1','r2','r3','k','kr','r2_dot','k_m_kr_dot','kr_inv'};
                for j=1:length(plotted)
                    ha(j)=subplot(3,3,j,'parent',hf); hold(ha(j),'on');
                    if j>6
                        xlabel(ha(j),'$t$ [$\mu$s]','interpreter','latex')
                    end
                    ylabel(ha(j),plotted{j},'interpreter','tex')
                    plot(ha(j),dim.st.time*1e6,dim.st.(plotted{j}),'k-','displayname','a priori');
                    plot(ha(j),dim.st.time(par.st.stable_time)*1e6,dim.st.(plotted{j})(par.st.stable_time),'rx','displayname','linear interpolation point');
                    plot(ha(j),dim.st.time(par.st.n_reconnection-par.st.rec_end_time)*1e6,dim.st.(plotted{j})(par.st.n_reconnection-par.st.rec_end_time),'bx','displayname','shift k-expression point');
                    h(j)=plot(ha(j),NaN,NaN,'m-','displayname','model output');
                end
                warning('off','MATLAB:linkaxes:RequireDataAxes')
                linkaxes(ha,'x');
            end
        end
        %% Switch in time
        if time>dim.st.time(end)
            % After sawtooth
            st.r1=0;
            st.r2=dim.st.r2(end);
            st.r3=st.r2;
            st.r1_dot=0;
            st.r2_dot=0;
            st.r3_dot=0;
            
            st.kr = 1;
            st.k = 0;
            
            st.k_m_kr_dot = 0;
        elseif time>par.st.t_reconnection
            % After reconnection (i.e. relaxation)
            st.r1=0;
            st.r2=dim.st.r2(end);
            st.r3=st.r2;
            st.r1_dot=0;
            st.r2_dot=0;
            st.r3_dot=0;
            
            st.kr=1;
            st.k=K(time);
            
            st.k_m_kr_dot = dim.st.k_m_kr_dot(end);
        elseif time<dim.st.time(par.st.stable_time)
            % Before the first stable point in reconnection
            st.r1       = dim.st.r1     (par.st.stable_time);
            st.r2       = dim.st.r2     (par.st.stable_time);
            st.r3       = dim.st.r3     (par.st.stable_time);
            st.r1_dot   = dim.st.r1_dot (par.st.stable_time);
            st.r2_dot   = dim.st.r2_dot (par.st.stable_time);
            st.r3_dot   = dim.st.r3_dot (par.st.stable_time);
            
            st.kr       = dim.st.kr     (par.st.stable_time);
            st.k        = dim.st.k      (par.st.stable_time);
            
            st.k_m_kr_dot = dim.st.k_m_kr_dot(par.st.stable_time);
        elseif time>dim.st.time(par.st.n_reconnection-par.st.rec_end_time)
            % Avoid problems with k at end of reconnection due to interpolation accuracy in r2 / r2 -dot
            st.r1=dim.st.r1_t(time);
            st.r2=R2(time);
            st.r3=sqrt(st.r2.^2-st.r1.^2);
            st.r1_dot=dim.st.r1_dot(1);
            st.r2_dot=R2_dot(time);
            st.r3_dot = (st.r2_dot.*st.r2-st.r1_dot.*st.r1)./st.r3;
            
            st.kr=(st.r2-st.r1)./(st.r2+st.r1);
            st.k=K(time);
            
            st.k_m_kr_dot = K_M_KR_dot(time);
        else
            % During reconnection (r1_dot = constant, k=kc)
            st.r1=dim.st.r1_t(time);
            st.r2=R2(time);
            st.r3=sqrt(st.r2.^2-st.r1.^2);
            st.r1_dot=dim.st.r1_dot(1);
            st.r2_dot=R2_dot(time);
            st.r3_dot = (st.r2_dot.*st.r2-st.r1_dot.*st.r1)./st.r3;
            
            st.kr=(st.r2-st.r1)./(st.r2+st.r1);
            st.k=-2*st.r1.*st.r2.* (st.r1_dot+st.r2_dot)./(st.r3.*st.r3_dot.*(st.r1+st.r2));
            
            st.k_m_kr_dot = K_M_KR_dot(time);
        end
        
        st.eta=st.r2-st.r1;
        st.kr_inv=1./st.kr;
        st.r3_sq=st.r3^2;
        
        
        if par.st.show_parameters && (...
            mod(time,1e-7)<=par.dt ||...
            (time>dim.st.time(end-5) && time < dim.st.time(end)+5e-6 &&   mod(time,1e-7)<=par.dt) ||...
            (abs(time-dim.st.time(par.st.n_reconnection-par.st.rec_end_time))<par.st.t_reconnection-dim.st.time(par.st.n_reconnection-par.st.rec_end_time))...
            )
            t_scaled=time*1e6;
            for j=1:length(plotted)
                X=get(h(j),'XData');
                Y=get(h(j),'YData');
                set(h(j),'XData',[X t_scaled],'YData',[Y st.(plotted{j})]);
            end

            xlim(ha(1),[max(0,t_scaled-5), t_scaled+10]);
            drawnow;
        end
    end
%% Sawtooth function 2 (space dependent parameters)
    function [Phi,Psi_star] = get_potential_flux(st,x,y,r_sq)
        %get_potential_flux [NESTED] Determines the Phi and Psi-functions
        %based on the sawtooth parameters for positions/matrixes x and y
        %and r^2!
        
        persistent PSI_MIN PSI_PLUS
        if isempty(PSI_MIN)
            PSI_MIN =griddedInterpolant(dim.st.r.^2,dim.st.psi_star,'linear');    % Flux in original core
            PSI_PLUS=griddedInterpolant(dim.st.r3(1:par.st.n_reconnection+1).^2,dim.st.psi_plus(1:par.st.n_reconnection+1),'linear');    % Flux in island
        end
        
        % Definition of untouched area
        expr_1 = r_sq    >=  st.r2.^2;       % Always untouched area
        
        % Pre-allocate Psi^* (always returned)
        Psi_star=zeros(size(x));
        
        % After sawtooth:
        if time > dim.st.time(end)
            % After sawtooth
            Phi=[];
            expr_3 = ~expr_1;
            Psi_star(expr_1) = PSI_MIN(r_sq(expr_1));
            Psi_star(expr_3) = PSI_PLUS(r_sq(expr_3));
            return
        end
        
        % Pre-allocate Phi
        Phi=zeros(size(x));
        
        if st.r1~=0 % During reconnection process
            zeta_sq = (x-st.eta).^2+y.^2;    % distance from center area 2 (original core)
            expr_2  =  zeta_sq <   st.r1.^2;       % Original core
            expr_3  = ~expr_1  &   ~expr_2;        % The rest, i.e. the island
            
            % General parameters
            R    =    sqrt((st.r2-x) .^2+y.^2);    % distance from x-point
            Theta_2=atan2(y(expr_2),st.r2-x(expr_2));
            Theta_3=atan2(y(expr_3),st.r2-x(expr_3));
            sin_Theta_3     = sin(Theta_3);
            cos_Theta_3     = cos(Theta_3);
            
            % r_plus calculation
            c=(st.kr_inv+st.k) - (R(expr_3)./(cos_Theta_3*st.r3)).^2;
            rho=c.*cos_Theta_3./(1 + sqrt(1+(st.k-st.kr)*c));
            r_plus_sq=(st.r3_sq*(rho.^2+sin_Theta_3.^2));
            
            % Psi^*  Helical flux
            Psi_star(expr_1) = PSI_MIN(r_sq(expr_1));
            Psi_star(expr_2) = PSI_MIN(zeta_sq(expr_2));
            Psi_star(expr_3) = PSI_PLUS(r_plus_sq);
            
            % Potential
            Phi(expr_2)=(st.r2_dot-st.r1_dot)*(R(expr_2).*sin(Theta_2));
            Phi(expr_3)= ...
                +0.5*Theta_3.*(r_plus_sq-st.r3_sq)*st.k_m_kr_dot...
                +sin_Theta_3.*(...
                st.r2_dot.*R(expr_3)...
                -(st.r1*st.r1_dot+st.r2*st.r2_dot)*cos_Theta_3...
                +st.r3_dot*st.r3.*rho);
            %+Theta_3.*(st.k-st.kc).*st.r3.*st.r3_dot... % This term should be zero
            
        else % Relaxation process
            expr_3 = ~expr_1; % The rest, i.e. the new core / island
            
            % General parameters
            R    =    sqrt((st.r2-x) .^2+y.^2);    % distance from x-point
            Theta_3=atan2(y(expr_3),st.r2-x(expr_3));
            sin_Theta_3   = sin(Theta_3);
            cos_Theta_3=    cos(Theta_3);
            
            % r_plus calculation
            c=(st.kr_inv+st.k) - (R(expr_3)./(cos_Theta_3*st.r3)).^2;
            rho=c.*cos_Theta_3./(1 + sqrt(1+(st.k-st.kr)*c));
            r_plus_sq=(st.r3_sq*(rho.^2+sin_Theta_3.^2));
            
            % Psi^*  Helical flux
            Psi_star(expr_1) = PSI_MIN(r_sq(expr_1));
            Psi_star(expr_3) = PSI_PLUS(r_plus_sq);
            
            % Potential
            Phi(expr_3)= ...
                +0.5*Theta_3.*(r_plus_sq-st.r3_sq)*st.k_m_kr_dot;
            %+sin_theta_3.*(...
            %st.r2_dot.*R(expr_3)...
            %-(st.r1*st.r1_dot+st.r2*st.r2_dot)*cos_theta_3...
            %+st.r3_dot*st.r3.*rho)...
            %+theta_3.*(st.k-st.kc).*st.r3.*st.r3_dot; % These terms should remain zero
        end
    end
end


%% B (3D)
function [BR,BZ,Bphi]     =B_3D_flux(x,X_ind,Z_ind)
global par dim maps
persistent F G H P LAST_SCHEME

% Find theta
theta = interpolate_theta_XZ(X_ind,Z_ind);  % Self-written linear interpolation inputs the indexes / slopes!

theta_ind=  theta               *dim.Dtheta_inv+1;
phi_ind=    mod(x(:,3,:),2*pi/dim.n3D.symm)  *dim.Dphi_inv+1;

%% 3D interpolation magnetic field

switch par.interp_scheme
    case 0 		%% SELF-WRITTEN LINEAR DIVERGENCE FREE / POTENTIAL BASED
        error('This interpolation scheme doesn''t support flux interpolation')
    case 1      %% SELF-WRITTEN LINEAR
        % Find 3D indexes based on theta,psi,phi
        psi_ind_3D  =interp2_XZ(X_ind,Z_ind,maps(1).psi_norm_3D_XZ,true(size(Z_ind.x)));
        [indexes,slopes] = interp_index_list_3D(dim.n3D.size_3D,theta_ind(:),psi_ind_3D,phi_ind(:));
        
        % Find fields
        BR      = lininterp3(maps(1).n3D.BR  ,indexes,slopes);
        BZ      = lininterp3(maps(1).n3D.BZ  ,indexes,slopes);
        Bphi    = lininterp3(maps(1).n3D.Bphi,indexes,slopes);
        
    case 2      %% INTERP LINEAR
        psi_ind_3D  =interp2(maps(1).psi_norm_3D_XZ,Z_ind,X_ind,'*linear');
        BR   = interp3(maps(1).n3D.BR  ,psi_ind_3D,theta_ind,phi_ind,'*linear');
        BZ   = interp3(maps(1).n3D.BZ  ,psi_ind_3D,theta_ind,phi_ind,'*linear');
        Bphi = interp3(maps(1).n3D.Bphi,psi_ind_3D,theta_ind,phi_ind,'*linear');
    case 3      %% BA_INTERP LINEAR
        psi_ind_3D  =ba_interp2(maps(1).psi_norm_3D_XZ,Z_ind,X_ind,'linear');
        BR    = ba_interp3(maps(1).B_3D,psi_ind_3D,theta_ind,phi_ind,'linear');      % Interpolate all fields at once
        BZ  =[];
        Bphi=[];
    case 4      %% GRIDDEDINTERPOLANT LINEAR
        if isempty(P)  || par.interp_scheme~=LAST_SCHEME
            P=griddedInterpolant(maps(1).psi_norm_3D_XZ,'linear');
            F=griddedInterpolant(maps(1).n3D.BR    ,'linear');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'linear');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'linear');
        end
        psi_ind_3D=P(X_ind,Z_ind);
        BR  =F(theta_ind,psi_ind_3D,phi_ind);
        BZ  =G(theta_ind,psi_ind_3D,phi_ind);
        Bphi=H(theta_ind,psi_ind_3D,phi_ind);
    case 5      %% INTERP CUBIC
        psi_ind_3D  =interp2(maps(1).psi_norm_3D_XZ,Z_ind,X_ind,'*cubic');
        BR   = interp3(maps(1).n3D.BR  ,psi_ind_3D,theta_ind,phi_ind,'*cubic');
        BZ   = interp3(maps(1).n3D.BZ  ,psi_ind_3D,theta_ind,phi_ind,'*cubic');
        Bphi = interp3(maps(1).n3D.Bphi,psi_ind_3D,theta_ind,phi_ind,'*cubic');
    case 6      %% BA_INTERP CUBIC
        psi_ind_3D  =ba_interp2(maps(1).psi_norm_3D_XZ,Z_ind,X_ind,'cubic');
        BR    = ba_interp3(maps(1).B_3D,psi_ind_3D,theta_ind,phi_ind,'cubic');      % Interpolate all fields at once
        BZ  =[];
        Bphi=[];
    case 7      %% GRIDDEDINTERPOLANT CUBIC
        if isempty(P) || par.interp_scheme~=LAST_SCHEME
            P=griddedInterpolant(maps(1).psi_norm_3D_XZ,'cubic');
            F=griddedInterpolant(maps(1).n3D.BR    ,'cubic');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'cubic');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'cubic');
        end
        psi_ind_3D=P(X_ind,Z_ind);
        BR  =F(theta_ind,psi_ind_3D,phi_ind);
        BZ  =G(theta_ind,psi_ind_3D,phi_ind);
        Bphi=H(theta_ind,psi_ind_3D,phi_ind);
    case 8      %% INTERP SPLINE
        psi_ind_3D  =interp2(maps(1).psi_norm_3D_XZ,Z_ind,X_ind,'*spline');
        BR   = interp3(maps(1).n3D.BR  ,psi_ind_3D,theta_ind,phi_ind,'*spline');
        BZ   = interp3(maps(1).n3D.BZ  ,psi_ind_3D,theta_ind,phi_ind,'*spline');
        Bphi = interp3(maps(1).n3D.Bphi,psi_ind_3D,theta_ind,phi_ind,'*spline');
    case 9      %% GRIDDEDINTERPOLANT SPLINE
        if isempty(P)  || par.interp_scheme~=LAST_SCHEME
            P=griddedInterpolant(maps(1).psi_norm_3D_XZ,'spline');
            F=griddedInterpolant(maps(1).n3D.BR    ,'spline');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'spline');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'spline');
        end
        psi_ind_3D=P(X_ind,Z_ind);
        BR  =F(theta_ind,psi_ind_3D,phi_ind);
        BZ  =G(theta_ind,psi_ind_3D,phi_ind);
        Bphi=H(theta_ind,psi_ind_3D,phi_ind);
end

%% Store last interpolation scheme
LAST_SCHEME=par.interp_scheme;

end

%% B total (3D)
function [BR,BZ,Bphi]     =B_3D_toroidal(x,X_ind,Z_ind)
%%B_total interpolation in R,Z,phi-grid
global par dim maps
persistent F G H LAST_SCHEME

% Indexes 2D
if nargin<3
    X_ind=((x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
    Z_ind=( x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
end
phi_ind=mod(x(:,3,:),2*pi/dim.n3D.symm)  *dim.Dphi_inv+1;

%% 3D interpolation magnetic field
switch par.interp_scheme
    case 0 		%% SELF-WRITTEN LINEAR DIVERGENCE FREE / POTENTIAL BASED
        R_inv=1./x(:,1,:);
        [BR,BZ,Bphi] = get_B_div_free_3D(X_ind,Z_ind,phi_ind,R_inv);
    case 1      %% SELF-WRITTEN LINEAR
        error('Self-written interpolation not available for R,Z,phi-grid')
    case 2      %% INTERP LINEAR
        BR   = interp3(maps(1).n3D.BR  ,Z_ind,X_ind,phi_ind,'*linear');
        BZ   = interp3(maps(1).n3D.BZ  ,Z_ind,X_ind,phi_ind,'*linear');
        Bphi = interp3(maps(1).n3D.Bphi,Z_ind,X_ind,phi_ind,'*linear');
    case 3      %% BA_INTERP LINEAR
        BR    = ba_interp3(maps(1).B_3D,Z_ind,X_ind,phi_ind,'linear');      % Interpolate all fields at once
        BZ  =[];
        Bphi=[];
    case 4      %% GRIDDEDINTERPOLANT LINEAR
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).n3D.BR    ,'linear');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'linear');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'linear');
        end
        BR  =F(X_ind,Z_ind,phi_ind);
        BZ  =G(X_ind,Z_ind,phi_ind);
        Bphi=H(X_ind,Z_ind,phi_ind);
    case 5      %% INTERP CUBIC
        BR   = interp3(maps(1).n3D.BR  ,Z_ind,X_ind,phi_ind,'*cubic');
        BZ   = interp3(maps(1).n3D.BZ  ,Z_ind,X_ind,phi_ind,'*cubic');
        Bphi = interp3(maps(1).n3D.Bphi,Z_ind,X_ind,phi_ind,'*cubic');
    case 6      %% BA_INTERP CUBIC
        BR    = ba_interp3(maps(1).B_3D,Z_ind,X_ind,phi_ind,'cubic');      % Interpolate all fields at once
        BZ  =[];
        Bphi=[];
    case 7      %% GRIDDEDINTERPOLANT CUBIC
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).n3D.BR    ,'cubic');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'cubic');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'cubic');
        end
        BR  =F(X_ind,Z_ind,phi_ind);
        BZ  =G(X_ind,Z_ind,phi_ind);
        Bphi=H(X_ind,Z_ind,phi_ind);
    case 8      %% INTERP SPLINE
        BR   = interp3(maps(1).n3D.BR  ,Z_ind,X_ind,phi_ind,'*spline');
        BZ   = interp3(maps(1).n3D.BZ  ,Z_ind,X_ind,phi_ind,'*spline');
        Bphi = interp3(maps(1).n3D.Bphi,Z_ind,X_ind,phi_ind,'*spline');
    case 9      %% GRIDDEDINTERPOLANT SPLINE
        if isempty(F)  || par.interp_scheme~=LAST_SCHEME
            F=griddedInterpolant(maps(1).n3D.BR    ,'spline');
            G=griddedInterpolant(maps(1).n3D.BZ    ,'spline');
            H=griddedInterpolant(maps(1).n3D.Bphi  ,'spline');
        end
        BR  =F(X_ind,Z_ind,phi_ind);
        BZ  =G(X_ind,Z_ind,phi_ind);
        Bphi=H(X_ind,Z_ind,phi_ind);
end

%% Store last interpolation scheme
LAST_SCHEME=par.interp_scheme;

    function [BR,BZ,Bphi] = get_B_div_free_3D(X_ind,Z_ind,phi_ind,R_inv)
        size_3D_intp0=dim.n3D.size_3D-1;
        
        % The index in the matrices / grid of the point 000 (i.e. where all
        % Beta / Gamma / Alpha's are stored
        X_ind_floor=floor(X_ind);
        Z_ind_floor=floor(Z_ind);
        phi_ind_floor=floor(phi_ind);
        
        % The local coordinates (from 0 to 1)
        r  =(X_ind  -X_ind_floor)  *dim.DX;
        z  =(Z_ind  -Z_ind_floor)  *dim.DZ;
        phi=(phi_ind-phi_ind_floor)*dim.Dphi;
        
        % Fast sub2ind
        index = X_ind_floor+(Z_ind_floor-1)*size_3D_intp0(1)+(phi_ind_floor-1)*size_3D_intp0(2)*size_3D_intp0(1);  %left  down  bottom corner
        index(isnan(index))=1;
        
        % Calculation steps
        Beta_r_r = maps(1).n3D.Alpha_r(index) + z.*maps(1).n3D.Gamma_z(index) + phi.*maps(1).n3D.Gamma_phi(index);
        Beta_z_z = maps(1).n3D.Alpha_z(index) - r.*maps(1).n3D.Gamma_r(index) - phi.*maps(1).n3D.Gamma_phi(index);
        
        BR   = (maps(1).n3D.Beta_r(index) + r .* Beta_r_r                     + z .* maps(1).n3D.Beta_r_z(index)    - phi .* maps(1).n3D.Beta_r_phi(index)).*R_inv;
        BZ   = (maps(1).n3D.Beta_z(index) + r .* maps(1).n3D.Beta_z_r(index)  + z .* Beta_z_z                       - phi .* maps(1).n3D.Beta_z_phi(index)).*R_inv;
        Bphi = maps(1).n3D.Beta_phi(index)+ r .* maps(1).n3D.Beta_phi_r(index)+ z .* (maps(1).n3D.Beta_phi_z(index) + r .* maps(1).n3D.Beta_phi_r_z(index)) - phi .* (Beta_z_z + Beta_r_r);
        
    end
end