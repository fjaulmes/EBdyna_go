function [ RMP_map,Bu,Bl,A] = BS_AUG_RMP_flux_surfaces(I_standard,RMP_map)
%BS_single_AUG_coil Determines the Magnetic field in R,Z,phi-cooirdinates
%   Detailed explanation goes here

delete(findall(0,'type','figure','tag','FLUX SURFACES'));
drawnow;

%% Load Coils
[ RMP_coils ] = BS_AUG_RMP_load_coils(1,1,0);
Bu=RMP_coils.Bu; Bl=RMP_coils.Bl; A=RMP_coils.A;
%% No input data
if nargin==0
    I_standard=4.8e3;
end

%% Definition of map dimensions
if exist('RMP_map','var')==0 || isempty(RMP_map)
    RMP_map.I=I_standard;
    [RMP_map,data_tokamak]=BS_AUG_flux_coordinates(129,129,15,RMP_map);
else
    if ~isfield(RMP_map,'I') || ~isfield(RMP_map,'R') || ~isfield(RMP_map,'Z') || ~isfield(RMP_map,'phi')
        error('Please provide RMP_map-struct contains all R,Z and phi-coordinates plus the current I for B-field calculation')
    end
end

%% Expressions
[~,expr_z0]=min(abs(RMP_map.theta(:,1,1)));
% [~,expr_z0]=min(abs(RMP_map.Z(1,1,:)));
[~,expr_y0]=min(abs(RMP_map.Y(:,1,1)));
expr_calc=RMP_map.R<(1.65+0.5) & RMP_map.R>(1.65-0.5);

%% Calculate B-field for each coil:
BS(1).Nfilament=1;
clc
disp('************************************************')
disp('Starting standard B-field calculation of each coil:')
disp('************************************************')

% Bu-coils:
for i=1:length(Bu)
    disp(['Bu coil ',num2str(i)])
    BS(i).Nfilament=1;  % Number of filaments of a coil
    BS(i).filament.Gamma = [Bu(i).X,Bu(i).Y,Bu(i).Z]'; % Coordinates coil
    BS(i).filament.I = RMP_map.I; % Current through coil
    BS(i).filament.dGamma = 1e9; % Maximum size of coil element dl (otherwise 2D interpolation of coil coordinates)
    
    RMP_map.Bu(i).BX=zeros(size(RMP_map.X));
    RMP_map.Bu(i).BY=zeros(size(RMP_map.X));
    RMP_map.Bu(i).BZ=zeros(size(RMP_map.X));
    % Calculate and store in RMP_map:
    RMP_map.Z=RMP_map.Z+0.07;   % Up the evaluation points by 0.07 since plasme is centered higher.
    [~,RMP_map.Bu(i).BX(expr_calc),RMP_map.Bu(i).BY(expr_calc),RMP_map.Bu(i).BZ(expr_calc)] = BS_calc_B(BS(i),RMP_map.X(expr_calc),RMP_map.Y(expr_calc),RMP_map.Z(expr_calc));
    RMP_map.Z=RMP_map.Z-0.07;   % Reduce Z express the vertical distance from center of plasma again
end
% Bl-coils:
for i=1:length(Bl)
    disp(['Bl coil ',num2str(i)])
    BS(i).Nfilament=1;  % Number of filaments of a coil
    BS(i).filament.Gamma = [Bl(i).X,Bl(i).Y,Bl(i).Z]'; % Coordinates coil
    BS(i).filament.I = RMP_map.I; % Current through coil
    BS(i).filament.dGamma = 1e9; % Maximum size of coil element dl (otherwise 2D interpolation of coil coordinates)
    
    RMP_map.Bl(i).BX=zeros(size(RMP_map.X));
    RMP_map.Bl(i).BY=zeros(size(RMP_map.X));
    RMP_map.Bl(i).BZ=zeros(size(RMP_map.X));
    % Calculate and store in RMP_map:
    RMP_map.Z=RMP_map.Z+0.07;   % Up the evaluation points by 0.07 since plasme is centered higher.
    [~,RMP_map.Bl(i).BX(expr_calc),RMP_map.Bl(i).BY(expr_calc),RMP_map.Bl(i).BZ(expr_calc)] = BS_calc_B(BS(i),RMP_map.X(expr_calc),RMP_map.Y(expr_calc),RMP_map.Z(expr_calc));
    RMP_map.Z=RMP_map.Z-0.07;   % Reduce Z express the vertical distance from center of plasma again
end
% A-coils:
for i=1:length(A)
    disp(['A coil ',num2str(i)])
    BS(i).Nfilament=1;  % Number of filaments of a coil
    BS(i).filament.Gamma = [A(i).X,A(i).Y,A(i).Z]'; % Coordinates coil
    BS(i).filament.I = RMP_map.I; % Current through coil
    BS(i).filament.dGamma = 1e9; % Maximum size of coil element dl (otherwise 2D interpolation of coil coordinates)
    
    RMP_map.A(i).BX=zeros(size(RMP_map.X));
    RMP_map.A(i).BY=zeros(size(RMP_map.X));
    RMP_map.A(i).BZ=zeros(size(RMP_map.X));
    % Calculate and store in RMP_map:
    RMP_map.Z=RMP_map.Z+0.07;   % Up the evaluation points by 0.07 since plasme is centered higher.
    [~,RMP_map.A(i).BX(expr_calc),RMP_map.A(i).BY(expr_calc),RMP_map.A(i).BZ(expr_calc)] = BS_calc_B(BS(i),RMP_map.X(expr_calc),RMP_map.Y(expr_calc),RMP_map.Z(expr_calc));
    RMP_map.Z=RMP_map.Z-0.07;   % Reduce Z express the vertical distance from center of plasma again
end

disp('************************************************')
disp('DONE!')
disp('************************************************')

%% Total field
BX.Bu=zeros(size(RMP_map.X)); BX.Bl=BX.Bu; BX.A=BX.Bu;
BY.Bu=zeros(size(RMP_map.X)); BY.Bl=BY.Bu; BY.A=BY.Bu;
BZ.Bu=zeros(size(RMP_map.X)); BZ.Bl=BZ.Bu; BZ.A=BZ.Bu;
for i=1:length(Bu)
    switch i
        case {1,2,5,6}
            Bu(i).sign=1;
        case {3,4,7,8}
            Bu(i).sign=-1;
    end
    % Add to field
    BX.Bu=BX.Bu+Bu(i).sign*RMP_map.Bu(i).BX;
    BY.Bu=BY.Bu+Bu(i).sign*RMP_map.Bu(i).BY;
    BZ.Bu=BZ.Bu+Bu(i).sign*RMP_map.Bu(i).BZ;
end

for i=1:length(Bl)
    switch i
        case {1,2,5,6}
            Bl(i).sign=1;
        case {3,4,7,8}
            Bl(i).sign=-1;
    end
    % Add to field
    BX.Bl=BX.Bl+Bl(i).sign*RMP_map.Bl(i).BX;
    BY.Bl=BY.Bl+Bl(i).sign*RMP_map.Bl(i).BY;
    BZ.Bl=BZ.Bl+Bl(i).sign*RMP_map.Bl(i).BZ;
end

for i=1:length(A)
    switch i
        case {1,2,5,6}
            A(i).sign=1;
        case {3,4,7,8}
            A(i).sign=-1;
    end
    % Add to field
    BX.A=BX.A+A(i).sign*RMP_map.A(i).BX;
    BY.A=BY.A+A(i).sign*RMP_map.A(i).BY;
    BZ.A=BZ.A+A(i).sign*RMP_map.A(i).BZ;
end

%% Calculation B-field
% In Cartesian
BX.all=BX.Bu+BX.Bl+BX.A;
BY.all=BY.Bu+BY.Bl+BY.A;
BZ.all=BZ.Bu+BZ.Bl+BZ.A;

% In cylindrical
BR=RMP_map.X./RMP_map.R.*BX.all+RMP_map.Y./RMP_map.R.*BY.all; % Component B in radial direction (to solenoid)
Bphi=-RMP_map.X./RMP_map.R.*BY.all+RMP_map.Y./RMP_map.R.*BX.all; % Component B in phi direction (to 

% In flux-coordinates

%% Integrate Bz for each R-coordinate
index_phi=2;
[RR, ZZ]=meshgrid(data_tokamak.scale_X+data_tokamak.R0,data_tokamak.scale_Z); %% R,Z-coordinates for contourplot

% Get Br,Bz,R,Z of calculation
RR2=RMP_map.R(:,:,index_phi);
ZZ2=RMP_map.Z(:,:,index_phi);
BBZZ=BZ.all(:,:,index_phi);
BBRR=BR(:,:,index_phi);

% Get Bz, Br for contourplot points with cubic interpolation
B_z_with=griddata(RR2(:),ZZ2(:),BBZZ(:),RR,ZZ,'cubic');
B_r=griddata(RR2(:),ZZ2(:),BBRR(:),RR,ZZ,'cubic');

% Field without RMP and with RMP  and total field
B_z_wo=data_tokamak.BZ_XZ';
B_z_tot=B_z_with+B_z_wo;
expr_outside_LCFS=isnan(B_z_with);

% Set to 0 for proper integration / flux calculation
B_z_wo(expr_outside_LCFS)=0;
B_z_with(expr_outside_LCFS)=0;
B_z_tot(expr_outside_LCFS)=0;

% Flux integration: psi=int(2*pi*r*B_z*dr)/2*pi
dFlux_with=B_z_tot.*RR;
dFlux_wo=B_z_wo.*RR;
Flux_Z_with=cumsum(dFlux_with(:,2:end).*diff(RR,1,2),2);
Flux_Z_with=[zeros(size(Flux_Z_with,1),1), Flux_Z_with];

Flux_Z_wo=cumsum(dFlux_wo(:,2:end).*diff(RR,1,2),2);
Flux_Z_wo=[zeros(size(Flux_Z_wo,1),1), Flux_Z_wo];

Flux_Z_with(expr_outside_LCFS)=NaN;
Flux_Z_wo(expr_outside_LCFS)=NaN;

%% Calculate nessecary B-maps
BX.all(~expr_calc)=NaN; BY.all(~expr_calc)=NaN; BZ.all(~expr_calc)=NaN;
X=permute(RMP_map.X,[3 2 1]); 	Y=permute(RMP_map.Y,[3 2 1]);	Z=permute(RMP_map.Z,[3 2 1]);

B_tot=sqrt(BX.all.^2+BY.all.^2+BZ.all.^2);
nabla_B=divergence(RMP_map.X,RMP_map.Y,RMP_map.Z,BX.all,BY.all,BZ.all); % Divergence

%% Contour of new flux surface
hf=figure('tag','FLUX SURFACES');
ha1=subplot(1,2,1,'parent',hf);
ha2=subplot(1,2,2,'parent',hf);
hold(ha1,'on'), grid(ha1,'on'), hold(ha2,'on'), grid(ha2,'on')


contour(RR(1,:),ZZ(:,1),Flux_Z_wo,100,'parent',ha1)
contour(RR(1,:),ZZ(:,1),Flux_Z_with,100,'parent',ha2)
axis(ha1,'equal','xy'),axis(ha2,'equal','xy'),  linkaxes([ha1 ha2],'xy')
colorbar('peer',ha1,'location','southoutside'), colorbar('peer',ha2,'location','southoutside')
xlabel(ha1,'$R$ [m]','interpreter','latex'), xlabel(ha2,'$R$ [m]','interpreter','latex')
ylabel(ha1,'$Z$ [m]','interpreter','latex'), ylabel(ha2,'$Z$ [m]','interpreter','latex')
title(ha1,['w/o RMP at $\phi=',num2str(RMP_map.phi(1,1,2)),'$-plane'],'interpreter','latex')
title(ha2,['with RMP at $\phi=',num2str(RMP_map.phi(1,1,2)),'$-plane'],'interpreter','latex')

R_LCFS=permute(RMP_map.R,[1 3 2]);
R_LCFS=R_LCFS(:,:,end);
Z_LCFS=permute(RMP_map.Z,[1 3 2]);
Z_LCFS=Z_LCFS(:,:,end);
plot(ha1,R_LCFS(:,2),Z_LCFS(:,2),'k','linewidth',3)
plot(ha2,R_LCFS(:,2),Z_LCFS(:,2),'k','linewidth',3)


%% Plot X,Y-plane figure (top view torus) (B)
hf=figure('tag','FLUX SURFACES');
ha1=axes('parent',hf);
hold(ha1,'on'), grid(ha1,'on')

B_tot_z0=permute(B_tot,[3 2 1]);
contourf(X(:,:,expr_z0),Y(:,:,expr_z0),abs(B_tot_z0(:,:,expr_z0)),'parent',ha1,'displayname','$\vec{B$')
c=get(gca,'children');
axis(ha1,'equal')
clear B_tot_z0
colorbar('peer',ha1)
xlabel(ha1,'$x$ [m]','interpreter','latex')
ylabel(ha1,'$y$ [m]','interpreter','latex')
title(ha1,'$B_z$ in [T]','interpreter','latex')

plot_coils(ha1,Bu,Bl,A)

view(ha1,3)
axis(ha1,'manual')
%% Poloidal cross section (movie)

            delete(c)

for i=1:size(RMP_map.phi,3)
    % Slice psi-theta-phi coordinates
    X2=RMP_map.psi(:,:,i);
    Y2=RMP_map.theta(:,:,i);
    Z2=RMP_map.phi(:,:,i);
    % Plot this
    hslicer = slice(RMP_map.psi,RMP_map.theta,RMP_map.phi,BZ.all,X2,Y2,Z2,'parent',ha1);
    set(hslicer,'EdgeColor','none')
    % Reset X-Y-Z-data
    set(hslicer,'XData',RMP_map.X(:,:,i),'YData',RMP_map.Y(:,:,i),'ZData',RMP_map.Z(:,:,i))
    colorbar('peer',ha1)
    
    drawnow;
    pause(0.3)
    if i~=2
        delete(hslicer)
    end
end
%% Adding steamlines on one flux surface
% General magnetic field values
%  Btot=sqrt(BX.all.^2+BY.all.^2+BZ.all.^2);
hf=figure('tag','FLUX SURFACES');
ha1=axes('parent',hf);

sr=R0+linspace(-0.3,0.3,5);
sz=linspace(-0.5,0.5,5);
[sr,sz]=meshgrid(sr,sz);
streamline(RR,ZZ,B_r,B_z_with,sr,sz,'parent',ha1);
axis(ha1,'equal','xy')
end
%% Function to create 2,2-subplot with certain name, X,Z-plane
function [ha1,ha2,ha3,ha4]=make_figure(name)
hf=figure('name',name,'tag','FLUX SURFACES');
ha1=subplot(2,2,1,'parent',hf);
ha2=subplot(2,2,2,'parent',hf);
ha3=subplot(2,2,3,'parent',hf);
ha4=subplot(2,2,4,'parent',hf);
hold(ha1,'on'), hold(ha2,'on'),hold(ha3,'on'),hold(ha4,'on')
for ha=[ha1 ha2 ha3 ha4]
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$z$ [m]','interpreter','latex')
end
end

%% Function to add coil shape in figure
function plot_coils(ha1,Bu,Bl,A)
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