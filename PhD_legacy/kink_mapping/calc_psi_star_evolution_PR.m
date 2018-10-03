% Wednesday April 25th 2012
% Evaluating the flux evolution through q=1 flux surface
% during a sawtooth collape
% according to initial and final profiles

load_tokamak_parameters;

Dr_avg=mean((dist_surf_PR_map(1:end-1,:)),1);

% Then find maximum value
xPsih_max=psi_rank_q1;
Psih=psi_star_initial_profile;
if psi_star_initial_profile(xPsih_max)>0
    psi_star_max=max(Psih);
    SIGN_PSIH=1;
else
    psi_star_max=min(Psih);
    Psih=-Psih;
    SIGN_PSIH=-1;
end

%Btor=Btor_initial_profile;
NRADIAL=Nradial;
volume_radial=volume_flux;
radial_r_value=radial_r_value_flux;
scale_X=(1:Nradial);



x_begin=round(0.5*(r_value_q1_mean/a)*NRADIAL);
% Find 0 value
xPsih_neg_pos=1;
xPsih_zero=interp1(Psih(psi_rank_q1:NRADIAL),psi_rank_q1:NRADIAL,0);
xMixing_radius=ceil(xPsih_zero);
% [epsilon xPsih_zero]=min(abs(Psih(psi_rank_q1:NRADIAL)));
% xPsih_zero=xPsih_zero+psi_rank_q1-1;

%Psih(xPsih_zero)=0;

Bstar_initial=zeros(1,NRADIAL);
Bstar_initial(1)=0;

for x=2:xMixing_radius
    %Dr_value=radial_r_value(x)-radial_r_value(x-1);
    Dr_value=Dr_avg(x);
    Bstar_initial(x)=2*((Psih(x)-Psih(x-1))/Dr_value)-Bstar_initial(x-1);
end




disp('Size of radial domain = (Position of mixing radius at x=)');
disp(xPsih_zero);
disp('Position of q=1 radius at x=');
disp(xPsih_max);

% disp('Position of mixing radius at r=');
% rPsih_zero=radial_r_value(xPsih_zero);
% disp(rPsih_zero);

%initial_Psimax_pos=radial_r_value(xPsih_max);



% reconstructing final psi curve is tricky
% we need to take less resolution in x

NPsih=xMixing_radius-xPsih_max+1;
Psih_values=zeros(1,NPsih);
volume_minus=zeros(1,NPsih);
volume_plus=zeros(1,NPsih);
volume_final=zeros(1,NPsih);

%define x_plus (more packed) values

for x=xMixing_radius:-1:xPsih_max
    psi_rank=(xMixing_radius-x)+1;
    Psih_values(psi_rank)=Psih(x);
    x_plus(psi_rank)=x;
    volume_plus(psi_rank)=volume_radial(x);
end
x_plus(1)=xMixing_radius;
x_plus(NPsih)=xPsih_max;

%initializing final Psih scaling to the initial one
%Psih_final_values=Psih_values;
xPsih_max=psi_rank_q1;
% if Psih(xPsih_max)>0
    Psih_values=min(Psih_values,Psih(xPsih_max));
% else
%     Psih_values=max(Psih_values,Psih(xPsih_max));
% end
Psih_values(NPsih-1)=0.5*(Psih_values(NPsih)+Psih_values(NPsih-2));

%find corresponding x_minus values
for n=2:NPsih-1
    x_minus(n)=interp1(Psih(1:xPsih_max),scale_X(1:xPsih_max),Psih_values(n));
    volume_minus(n)=interp1(Psih(1:xPsih_max),volume_radial(1:xPsih_max),Psih_values(n));
end
x_minus(1)=1;
x_minus(NPsih)=xPsih_max;
volume_minus(1)=volume_radial(1);
volume_minus(NPsih)=volume_radial(NPsih);



%building r and x values for final Psih function
x_psi_final(1)=xMixing_radius;
for n=2:NPsih-1
    volume_final(n)=volume_plus(n)-volume_minus(n);
    x_psi_final(n)=interp1(volume_radial(1:xMixing_radius),scale_X(1:xMixing_radius),volume_final(n));
end
x_psi_final(NPsih)=1;


    
% building values for final psi function in increasing order
for y=1:NPsih
    x_final_interp(y)=x_psi_final(NPsih-y+1);
    Psih_final_interp(y)=Psih_values(NPsih-y+1);
end

% And now we find back the poloidal field after crash
% To do this, we interpolate linearly between the points
Psih_final=interp1(x_final_interp,Psih_final_interp,(1:xMixing_radius),'cubic');


for x=xMixing_radius+1:NRADIAL
    %q_final(x)=q_initial(x);
    Psih_final(x)=Psih(x);
end




Bstar_final=zeros(1,NRADIAL);
Bstar_final(1)=0;

Bstar_PR_final=Bstar_PR_map;

for p=1:NP
    Bstar_PR_final(p,1)=0;
    Bstar_PR_final(p,2)=(Psih_final(3)-Psih_final(1))/(dist_surf_PR_map(p,3)+dist_surf_PR_map(p,2));
    
    for x=3:xMixing_radius
        Bstar_PR_final(p,x)=(Psih_final(x+1)-Psih_final(x-1))/(dist_surf_PR_map(p,x+1)+dist_surf_PR_map(p,x));
    end
end
Bstar_PR_final=Bstar_PR_final./Rpos_PR_map;

Bstar_final=mean(Bstar_PR_final(1:NP-1,:),1);

%Bstar_final(2)=(Psih_final(3)-Psih_final(1))/(2*DX);
Bstar_final(NRADIAL-1)=Bstar_final(NRADIAL-2);
Bstar_final(NRADIAL)=Bstar_final(NRADIAL-1);

Bpol_final=BHpol_profile+Bstar_final;
q_final(1:xMixing_radius)=BHpol_profile(1:xMixing_radius)./Bpol_final(1:xMixing_radius);
q_final(xMixing_radius+1:NRADIAL)=q_initial(xMixing_radius+1:NRADIAL);



% Total magnetic field values
for x=1:NRADIAL
    Btot_final(x)=sqrt(Bpol_final(x)^2+Btor(x)^2);
    Btot_initial(x)=sqrt(Bpol_initial(x)^2+Btor(x)^2);
end




save_kadomtsev_profiles;