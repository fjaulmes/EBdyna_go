function [w ]=deviate_velocity_vector(v,v_norm,delta_omega)
%%time_step_integration_GT_eq_struct
% Multifunctional timestep with BORIS or Fabien method.
N_calc=length(v_norm);

% deviation=v_norm*0;

vu=bsxfun(@times,1./v_norm,v);
uu=zeros(N_calc,3);
uu(:,1)=sqrt(1./(1+(vu(:,1)./vu(:,3)).^2));  % x-component
uu(:,3)=-(vu(:,1)./vu(:,3)).*uu(:,1);        % phi-component

% Other perpendicular vector by cross product
wu=cross(uu,vu,2);
eR_conv=[vu(:,1) uu(:,1) wu(:,1)];
eZ_conv=[vu(:,2) uu(:,2) wu(:,2)];
ephi_conv=[vu(:,3) uu(:,3) wu(:,3)];

% new vector is w
w_norm=v_norm;

% closed intervall on 0 : isotropic diffusion
dispersion_angle=(randi([0, 2^52-1],N_calc,1) / 2^52)*2*pi;


delta_norm=v_norm.*delta_omega;
% in the plane orthogonal to v
delta_Y=delta_norm.*sin(dispersion_angle);
delta_Z=delta_norm.*cos(dispersion_angle);

% small deviation : coordinates in frame of current velocity
wX=sqrt(w_norm.^2-(delta_norm).^2);
wY=delta_Y;
wZ=delta_Z;

w_vec=[wX wY wZ];

wR=sum(w_vec.*eR_conv,2);
wZ=sum(w_vec.*eZ_conv,2);
wphi=sum(w_vec.*ephi_conv,2);

w=[wR wZ wphi];
% w_v=[wR-v(:,1) wZ-v(:,2) wphi-v(:,3)]; % required for calculation precision!
% 
% deviation = acos(1+(sum(w_v.*v,2))./(v_norm.*w_norm));

% v=w;
return