function [ n3D_map,Z_correction ] = BS_AUG_get_grid( n3D_map )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global par
narginchk(1,1)
Z_correction=0; % Set to 0 initially

switch par.example
    case 0
        %% NORMAL ROUTINE
        switch par.coord_sys
            case 'flux'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Make sure these psi and theta are 2^L+1 (L positive integer)
                nr_psi  =513;
                nr_theta=513;
                nr_phi  =1025;
                
                %% Load data and define size                
                par.size_total=[nr_theta nr_psi nr_phi];
                par.indexes=get_expr_job(par.size_total,par.N_PROCESS,par.PROCESS_NUMBER,'ind'); %Parallize on phi
                dim=load('../data_tokamak/motions_map_dimensions.mat','Z_axis');
                
                [n3D_map,data_tokamak]=BS_AUG_flux_coordinates(n3D_map); % Determine with flux coordinates (psi,theta,phi)
                
                Z_correction=double(data_tokamak.z_axis)-dim.Z_axis;
                %% Determine indexes used for finesse mesh
                data_tokamak.theta_index=interp1(1:data_tokamak.size(1),linspace(1,data_tokamak.size(1),nr_theta),'nearest');
                data_tokamak.psi_index=interp1(1:data_tokamak.size(2),linspace(1,data_tokamak.size(2),nr_psi),'nearest');
                if any(diff(data_tokamak.theta_index)==0) || any(diff(data_tokamak.psi_index)==0)
                    error('Flux coordinate grid not properly chosen')
                end
            case 'toroidal'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dim=load('../data_tokamak/motions_map_dimensions.mat', 'scale_X', 'scale_Z','Z_axis');
                
                nr_X=1*numel(dim.scale_X);
                nr_Z=1*numel(dim.scale_Z);
                nr_phi=513;
                
                par.size_total=[nr_X nr_Z nr_phi];
                par.indexes=get_expr_job(par.size_total,par.N_PROCESS,par.PROCESS_NUMBER,'ind'); %Parallize on phi
                                
                % Load FINESSE equilibrium data and make coordinate system
                [n3D_map,data_tokamak]=BS_AUG_toroidal_coordinates(n3D_map);    % Determine with flux coordinates (psi,theta,phi)
                Z_correction=double(data_tokamak.z_axis)-dim.Z_axis;                       % force the Z-offset to double precision. Then interpolation in grid has also double precision
                % NOTE: distance geometrical center to axis : z_axis
                %       distance Z=0 to axis                : Z_axis
                %    -> distance geometrical center to Z=0  : z_axis - Z_axis
            otherwise
                error('Coordinate system undefined')
        end
        %% Extra grid if vectorpotential
        if par.determine_vector_potential
            % Positions in extra positions for in phi-coordinate for dA/dphi-terms to provide
            % canonical toroidal angular momentom            
            avg_angle=2*pi/(nr_phi-1);          % difference in toroidal angle
            par.delta_phi2=avg_angle/20;        % Angle difference around point in grid 1
            phi_arounds=par.delta_phi2*[-2 -1 1 2]; % Angles around point 1
            n3D_map.R2=repmat(n3D_map.R,1,4);
            n3D_map.Z2=repmat(n3D_map.Z,1,4);
            phi2=bsxfun(@plus,phi_arounds,repmat(n3D_map.phi,1,4));
            
            % Cartesian coordinates:
            n3D_map.X2= n3D_map.R2.*cos(phi2);
            n3D_map.Y2=-n3D_map.R2.*sin(phi2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Examples (remake a grid)
    case 1
        % Cartesian cross section at y=0
        par.determine_vector_potential=false;
        
        x=linspace(1.65-0.5,1.65+0.5,1e2);
        y=[0 eps];
        z=linspace(-0.85,0.85,1e1);
        par.expr_y0=y==0;
        
        [n3D_map.X,n3D_map.Y,n3D_map.Z]=meshgrid(x,y,z);
        n3D_map.R=sqrt(n3D_map.X.^2+n3D_map.Y.^2);
        n3D_map.phi=mod(atan(-n3D_map.Y./n3D_map.X)+(sign(n3D_map.X)-1)*pi/2,2*pi);
        par.size_total=size(n3D_map.X);
        
        expr_job=get_expr_job(par.size_total,par.N_PROCESS,par.PROCESS_NUMBER); %Parallize
        expr_calc=n3D_map.R<(1.65+0.5) & n3D_map.R>(1.65-0.5);
        
        expr_job=expr_job & expr_calc;
        par.indexes=find(expr_job);
    case 2
        % Toroidal cross section at z=0
        par.determine_vector_potential=false;
        
        z=unique([0 linspace(-0.8,0.87,2)]);
        R=linspace(1.65-0.5,1.65+0.5,1e2);
        phi=linspace(0,2*pi,1e2);
        [n3D_map.R,n3D_map.phi,n3D_map.Z]=meshgrid(R,phi,z);
        n3D_map.X=n3D_map.R.*cos(n3D_map.phi);
        n3D_map.Y=-n3D_map.R.*sin(n3D_map.phi);
        
       par.size_total=size(n3D_map.X);
        
        expr_job=get_expr_job(par.size_total,par.N_PROCESS,par.PROCESS_NUMBER); %Parallize
        par.indexes=find(expr_job);
    case 3
        par.determine_vector_potential=false;
        % Cartesian cross section in poloidal plane
        z=unique([0 linspace(-0.8,0.87,2)]);
        
        x=linspace(-1.65-0.5,1.65+0.5,1e2);
        y=x;
        [n3D_map.X,n3D_map.Y,n3D_map.Z]=meshgrid(x,y,z);
        n3D_map.R=sqrt(n3D_map.X.^2+n3D_map.Y.^2);
        n3D_map.phi=mod(atan(-n3D_map.Y./n3D_map.X)+(sign(n3D_map.X)-1)*pi/2,2*pi);
        
        par.size_total=size(n3D_map.X);
        expr_job=get_expr_job(par.size_total,par.N_PROCESS,par.PROCESS_NUMBER); %Parallize
        par.indexes=find(expr_job);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reduce RMP_map to job size
if par.example~=0
    n3D_map=Reduce_RMP_map_job(par.size_total,n3D_map,expr_job);
    if par.determine_vector_potential
        n3D_map=Reduce_RMP_map_job(size_total2,n3D_map,expr_job_2);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fuction to reduce size of RMP (NOTE, THEY TURN INTO VECTORS)
function struct=Reduce_RMP_map_job(size_map,struct,expr_job)
fnames=fieldnames(struct);
for i=1:length(fnames)
    if isequal(size_map,size(struct.(fnames{i})))
        struct.(fnames{i})=struct.(fnames{i})(expr_job);
    end
end
end
