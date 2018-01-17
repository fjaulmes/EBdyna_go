
FMAX=size(psi_star_2D_evol_lin,1);

f_prev=max(f-2,1);
f0=f_prev;
initialize_rotated_PR_maps;
psi_star_PR_map_prev_prev=psi_star_PR_map;

f_next=min(f+2,FMAX);
f0=f_next;
initialize_rotated_PR_maps;
psi_star_PR_map_next_next=psi_star_PR_map;


f_prev=max(f-1,1);
f0=f_prev;
initialize_rotated_PR_maps;
psi_star_PR_map_prev=psi_star_PR_map;
if (f-1)>0
  time_evol_prev=time_scale_lin(f_prev);
else
  time_evol_prev=-time_scale_lin(f+1);
end
  
f_next=min(f+1,FMAX);
f0=f_next;

initialize_rotated_PR_maps;
psi_star_PR_map_next=psi_star_PR_map;

if (f+1)<size(time_scale_lin,2)
  time_evol_next=time_scale_lin(f_next);
else
  time_evol_next=2*time_scale_lin(f)-time_scale_lin(f-1);
end

% Here should be uniform time steps !!

COLLAPSE_DURATION=4e-4;
DT=COLLAPSE_DURATION*0.5*(time_evol_next-time_evol_prev);

psi_star_dot_PR_map=(1/12)*(-psi_star_PR_map_next_next+psi_star_PR_map_prev_prev)+(2/3)*(psi_star_PR_map_next-psi_star_PR_map_prev);
psi_star_dot_PR_map=psi_star_dot_PR_map/DT;


% psi_star_dot_RZ_map=griddata(finesse_data_X,finesse_data_Z,psi_star_dot_data,XX,ZZ,'cubic');
% psi_star_dot_RZ_map=psi_star_dot_RZ_map';
% psi_star_dot_RZ_map(isnan(psi_star_dot_RZ_map))=0;



%disp('psi_star_dot_RZ_map done ...')