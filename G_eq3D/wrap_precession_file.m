


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
%PART_POP=find(~alphas_ejected);
%PART_POP=find((~alphas_ejected).*(~isnan(vpll_output(end,:)))');
PART_POP=find((~alphas_ejected).*(~isnan(alphas_vpll)));

TIME_STAMP_PRECISION=1
TINI=TIME_STAMP_PRECISION;
NB_TIME_STAMPS=size(posX_output,1)
NB_TIME_STAMPS_NEW=round(1.0*NB_TIME_STAMPS)
time_scale_new=time_scale(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW);
time_scale=time_scale_new;


% pos_x=posX_output(~alphas_ejected);
% pos_z=posZ_output(~alphas_ejected);

X_output_pop=posX_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,PART_POP);
clear posX_output 
Z_output_pop=posZ_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,PART_POP);
clear posZ_output 
phi_output_pop=phipos_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,PART_POP);
clear phipos_output
vpll_output_pop=vpll_output(TINI:TIME_STAMP_PRECISION:NB_TIME_STAMPS_NEW,PART_POP);
clear vpll_output


% save('extra_precession_data.mat','INPUTNAME','mHe','ZHe','density_part_ratio','drX_map','drZ_map','X_output_pop','Z_output_pop','phi_output_pop','vpll_output_pop','alphas_Ekin','alphas_mm','alphas_pphi0','time_scale');
