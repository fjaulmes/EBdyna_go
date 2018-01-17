%------------------------------------------------------------
% Transform maps in PR coordinates so that they 
% have a rotation matching their toroidal position
%------------------------------------------------------------


psi_star_PR_map=zeros(NP,Nradial);
psi_star_PR_map_raw=zeros(NP,Nradial);

for (p=1:NP)
    psi_star_PR_map_raw(p,:)=psi_star_initial.*ones(1,Nradial);
end

%phi=0;

NP_half=(NP+1)/2;
psi_star_2D=zeros(size_r,NP_half);
psi_star_2D=squeeze(psi_star_2D_evol_lin(f0,:,:));
%for (p=1:NP_half)
%    for (r=1:size_r)
%        psi_star_2D(r,p)=psi_star_2D_evol_lin(f0,r,p);
%    end
%end

psi_star_PR_map_raw(1:NP_half,1:size_r)=psi_star_2D';
psi_star_PR_map_raw(NP:-1:NP_half,1:size_r)=psi_star_2D';


NP_begin=round(abs(phi)/Domega)+1;
NP_middle_low=NP_half-NP_begin+1;
NP_half=(NP+1)/2;
NP_middle_high=NP_middle_low+NP_half-1;

% psi_star_PR_map(1:NP_middle_low,1:size_r)=psi_star_2D(:,NP_begin:NP_half)';
% psi_star_PR_map(NP_middle_low:NP_middle_high,1:size_r)=psi_star_2D(:,NP_half:-1:1)';
% psi_star_PR_map(NP_middle_high:NP,1:size_r)=psi_star_2D(:,1:NP_begin)';

if (phi<=0)

    psi_star_PR_map(1:NP_middle_high,:)=psi_star_PR_map_raw(NP_begin:NP,:);
    psi_star_PR_map(NP_middle_high:NP,:)=psi_star_PR_map_raw(1:NP_begin,:);
    
%     Bstar_PR_map_rotated(1:NP_middle_high,:)=Bstar_PR_map(NP_begin:NP,:);
%     Bstar_PR_map_rotated(NP_middle_high:NP,:)=Bstar_PR_map(1:NP_begin,:);
%     
    NP_initial=NP_middle_low;
    NP_middle=NP_middle_high;
else

    psi_star_PR_map(1:NP_begin,:)=psi_star_PR_map_raw(NP_middle_high:NP,:);
    psi_star_PR_map(NP_begin:NP,:)=psi_star_PR_map_raw(1:NP_middle_high,:);
    
%     Bstar_PR_map_rotated(1:NP_begin,:)=Bstar_PR_map(NP_middle_high:NP,:);
%     Bstar_PR_map_rotated(NP_begin:NP,:)=Bstar_PR_map(1:NP_middle_high,:);
% 
    NP_initial=NP_begin;
    NP_middle=NP_middle_low;
end

psi_star_PR_map_raw=psi_star_PR_map;



