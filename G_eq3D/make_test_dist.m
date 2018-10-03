function [ x,vpll,vperp ] = make_test_dist( const,dim,type,mass )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

switch type
    case 1
        Nalphas_simulated=128;
        Nalphas_half=round(0.5*Nalphas_simulated)-1;
        Nalphas_fourth=round(0.25*Nalphas_simulated)-1;
        alphas_pos_z=zeros(1,Nalphas_simulated);
        
        alphas_pos_z(1:Nalphas_fourth+1)=0.3*((0:Nalphas_fourth)/Nalphas_fourth);
        alphas_pos_z(Nalphas_fourth+2:Nalphas_half+1)=-alphas_pos_z(1:Nalphas_fourth+1)-0.3;
        alphas_pos_z(1:Nalphas_fourth+1)=alphas_pos_z(1:Nalphas_fourth+1)+0.3;
        alphas_pos_z=alphas_pos_z+0.001;
        
        alphas_pos_x=alphas_pos_z*0;
        alphas_pos_x(round(0.5*Nalphas_simulated)+1:Nalphas_simulated)=0.9*((0:Nalphas_half)/Nalphas_half-0.5)';
        alphas_pos_x=alphas_pos_x-0.001;
        alphas_pos_z=alphas_pos_z' + dim.Z_axis;
        alphas_pos_x=alphas_pos_x';
        
        alphas_pos_phi=0*alphas_pos_x;
        
        alphas_Ekin=(90e3)*ones(Nalphas_simulated,1);
        alphas_vtot=sqrt(alphas_Ekin*2*const.eV/mass);
        vpll=(alphas_vtot*0.2);
        vpll(1:2:Nalphas_simulated)=-vpll(1:2:Nalphas_simulated);
 
        alphas_Epll=0.5*(mass/const.eV)*vpll.^2;
        alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
      
        %% Storing initial values and pre-allocating

        vperp=sqrt(2*alphas_Eperp*const.eV/mass); %Perpendicular velocity
        
        x=[alphas_pos_x+dim.R0,alphas_pos_z,alphas_pos_phi];
    case {2,3}
        %% Distribution homogenously distributed in ST-region with low energy
        N_particles=128;
       
        mes_size=ceil(sqrt(N_particles));
        N_particles=mes_size^2;
        Ekin=0.04 * ones(N_particles,1);
         
        x=linspace(1.50,1.85,mes_size);
        z=linspace(-0.22,0.22,mes_size);
        [x,z]=ndgrid(x,z);
        phi=zeros(size(x));
        
        x=cat(2,x(:),z(:),phi(:));
        
        % with E=Eperp of E=Epll
        % case 2: E=Eperp -> E times B drift
        % case 3: E=Epll  -> psi^* contours when E=0
        vtot=sqrt(Ekin*2*const.eV/mass);
        vperp=vtot;
        vpll=zeros(size(vtot));  
        if type==3
            vpll=vtot; vperp(:)=0;
        end
    case 4
        %% Test with few particles in deep core, should move with r1_dot - r2_dot
        N_particles=128;
        x=zeros(N_particles,3);
        x(:,1)=dim.R0+dim.X_axis;  % Put all on magnetic axis
        x(2:end,1:2)=x(2:end,1:2)+rand(N_particles-1,2)*1e-3; % Add small deviation
        
        Ekin=0.04*ones(N_particles,1);
        vtot=sqrt(Ekin*2*const.eV/mass);
        vperp=vtot;
        vpll=zeros(size(vperp));          
    case 5
        %% Distribution homogenously distributed in ST-region with low energy
        N_particles=128;
       
        mes_size=ceil(sqrt(N_particles));
        N_particles=mes_size^2;
        Ekin=2.8e3*ones(N_particles,1);
         
        x=linspace(1.50,1.85,mes_size);
        z=linspace(-0.22,0.22,mes_size);
        [x,z]=ndgrid(x,z);
        phi=zeros(size(x));
        
        x=cat(2,x(:),z(:),phi(:));
        
        vtot=sqrt(Ekin*2*const.eV/mass);
        rng(1); % seed
        prll_chance=1-2*rand(N_particles,1);
        vpll=vtot.*prll_chance;
        vperp=sqrt(vtot.^2-vpll.^2);
    case 9
        maps2=load('../data_tokamak/flux_geometry.mat');
        
        N_particles = 3;
        Ekin=100e3*ones(N_particles,1);
        vtot=sqrt(Ekin*2*const.eV/mass);
                
        psi_overline=1-dim.psi_scale./max(dim.psi_scale);
        psi_ind=interp1(psi_overline,1:dim.NB_PSI,0.2);
        
%         pphi=-const.eV*interp1(1:dim.NB_PSI,dim.psi_scale,psi_ind);
%         
%         R_scale=dim.R0+interp2(maps2.X_PR_map,1:513,1);
%         vpll_scale=(pphi+const.eV*dim.psi_scale)./(mass*R_scale);
%         x=interp1(vpll_scale,R_scale,vtot.*[-1 0 1]');

        
        theta_ind=[385 449 385]';
        
        x=dim.R0+interp2(maps2.X_PR_map,psi_ind,theta_ind);
        z=interp2(maps2.Z_PR_map,psi_ind,theta_ind);
        
        vpll=vtot.*[-0.7 0 0.7]';
        vperp=sqrt(vtot.^2-vpll.^2);
        phi=zeros(size(x));
        
        x=cat(2,x(:),z(:),phi(:));
        
       
end
end