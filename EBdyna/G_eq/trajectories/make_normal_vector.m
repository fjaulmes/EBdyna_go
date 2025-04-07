    function [N,b,gyro_angle]=make_normal_vector(input,x,N_simulated)
%% Find normal vector b
% Determine direction B with inclusion of RMP-field
if size(x,1)<2; error('Code errors likely for less than 2 particles per job. Make sure there are more particles than jobs'); end;
[~,B]=B_interpolation(x,'2D');
Bfield_inv=dot(B,B,2).^-0.5; % inverse magnitude
b=bsxfun(@times,Bfield_inv,B);

%% Make 2 perpendicular vectors
% Perpendicular unit vector with perpendicular to b and eZ (unit in Z-direction)
u=zeros(N_simulated,3);
u(:,1)=sqrt(1./(1+(b(:,1)./b(:,3)).^2)); % x-component
u(:,3)=-(b(:,1)./b(:,3)).*u(:,1);        % phi-component

% Other perpendicular vector by cross product
w=cross(u,b,2);

%% Make vector normal vector of perpendicular velocity
% normal vector N = cos(phase)* u + sin(phase)* w
seed=1;
rng(seed);
gyro_angle=rand(size(x,1),1)*2*pi;
% gyro_angle=rand(2000,1)*2*pi;
% sub_ind=mod(input.particle_nr-1,2000)+1;
% gyro_angle=gyro_angle(sub_ind);

N=   bsxfun(@times,cos(gyro_angle),u)...
    +bsxfun(@times,sin(gyro_angle),w);
end

