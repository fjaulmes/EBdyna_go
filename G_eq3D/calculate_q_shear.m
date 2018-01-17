load('../data_tokamak/q_profile.mat');
load('../data_tokamak/flux_geometry.mat');

d_radial=0.01;
radial_scale=0:d_radial:1;
q_profile_radial=interp1(radial_r_value_flux,q_initial_profile,radial_scale);


for(r=2:100)
    shear_profile(r)=(q_profile_radial(r+1)-q_profile_radial(r-1))/(2*d_radial);
end

shear_profile(1)=shear_profile(2);
shear_profile(101)=shear_profile(100);
% shear_profile=shear_profile';

shear_profile=(radial_scale./q_profile_radial).*shear_profile;

plot(radial_scale,shear_profile)

[val_min q1_rank ]=min(abs(q_profile_radial-1));
disp('shear_profile(q1_rank)=');
disp(shear_profile(q1_rank))

shear_profile_inital=interp1(radial_scale,shear_profile,radial_r_value_flux);
