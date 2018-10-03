


finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
finesse_mesh=[finesse_data_X  finesse_data_Z];
[finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
finesse_mesh_dtri=DelaunayTri(finesse_mesh);
dr_data=reshape(dr_X_PR_map(:,1:size_r),NP*size_r,1);
drX_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,'quadratic');
drX_map=drX_map';
dr_data=reshape(dr_Z_PR_map(:,1:size_r),NP*size_r,1);
drZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,'quadratic');
drZ_map=drZ_map';


load(INPUTNAME);
disp('NUMBER OF EJECTED PARTICLES')
disp(length(find(alphas_ejected)));

posX_output=Xpos_gc_output;
posZ_output=Zpos_gc_output;
% EJECTED_POP=find(alphas_ejected);

TIME_STAMP_PRECISION=1
TINI=TIME_STAMP_PRECISION;
NB_TIME_STAMPS=size(posX_output,1)
NB_TIME_STAMPS_NEW=round(1.0*NB_TIME_STAMPS)
time_scale_new=time_scale(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW);
time_scale=time_scale_new;

X_output_pop=posX_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,:);
clear posX_output
Z_output_pop=posZ_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,:);
clear posZ_output
phi_output_pop=phipos_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,:);
clear phipos_output
vpll_output_pop=vpll_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,:);
clear vpll_output

% this is done at the end of GT_many_precession_evolution_eq
%alphas_Ekin=alphas_Ekin(PART_POP);
%alphas_mm=alphas_mm(PART_POP);


disp('SAVING ALL PARTICLES (including post crash eq losses):')
length(alphas_Ekin)
length(vpll_output_pop)

save('extra_precession_data.mat','INPUTNAME','FILENUMBER','mHe','ZHe','density_part_ratio','drX_map','drZ_map','alphas_ejected','X_output_pop','Z_output_pop','phi_output_pop','vpll_output_pop','alphas_Ekin','alphas_mm','alphas_pphi0','time_scale');
