frame_rank_precise=interp1(time_scale_precise*tau_cr/FAST_SAWTOOTH,1:1001,time);
frame_rank_prev=10*floor((frame_rank_precise-1)/10)+1;
frame_rank_next=10*ceil((frame_rank_precise-1)/10)+1;

if frame_rank_next~=frame_rank_next_prev
    

    frame_rank=frame_rank_prev;
    vExB_X_map_phi_prev=vExB_X_map_phi_next;
    vExB_Z_map_phi_prev=vExB_Z_map_phi_next;
    vExB_phi_map_phi_prev=vExB_phi_map_phi_next;
    vD_X_map_phi_prev=vD_X_map_phi_next;
    vD_Z_map_phi_prev=vD_Z_map_phi_next;
    vD_phi_map_phi_prev=vD_phi_map_phi_next;
    bX_map_phi_prev=bX_map_phi_next;
    bZ_map_phi_prev=bZ_map_phi_next;
    bphi_map_phi_prev=bphi_map_phi_next;
    
    frame_rank=frame_rank_next;
    run('load_vExBmap_frame_rank');
    run('load_vD_b_frame_rank');
    vExB_X_map_phi_next=FAST_SAWTOOTH*vExB_X_map_phi;
    vExB_Z_map_phi_next=FAST_SAWTOOTH*vExB_Z_map_phi;
    vExB_phi_map_phi_next=FAST_SAWTOOTH*vExB_phi_map_phi;
    vD_X_map_phi_next=vD_X_map_phi;
    vD_Z_map_phi_next=vD_Z_map_phi;
    vD_phi_map_phi_next=vD_phi_map_phi;
    bX_map_phi_next=bX_map_phi;
    bZ_map_phi_next=bZ_map_phi;
    bphi_map_phi_next=bphi_map_phi;
    
    disp('frame_rank_precise = ');
    disp(frame_rank_precise);
    
end

vExB_X_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vExB_X_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vExB_X_map_phi_next;
vExB_Z_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vExB_Z_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vExB_Z_map_phi_next;
vExB_phi_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vExB_phi_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vExB_phi_map_phi_next;

vD_X_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vD_X_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vD_X_map_phi_next;
vD_Z_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vD_Z_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vD_Z_map_phi_next;
vD_phi_map_phi=0.1*(frame_rank_next-frame_rank_precise)*vD_phi_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*vD_phi_map_phi_next;

bX_map_phi=0.1*(frame_rank_next-frame_rank_precise)*bX_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*bX_map_phi_next;
bZ_map_phi=0.1*(frame_rank_next-frame_rank_precise)*bZ_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*bZ_map_phi_next;
bphi_map_phi=0.1*(frame_rank_next-frame_rank_precise)*bphi_map_phi_prev+0.1*(frame_rank_precise-frame_rank_prev)*bphi_map_phi_next;




%frame_rank_prev_prev=frame_rank_prev;
frame_rank_next_prev=frame_rank_next;


% vExB_X_PR_map=zeros(NP,size_r);
% vExB_X_XZ_map=zeros(size_X,size_Z);
% vExB_X_XZ_map_phi=zeros(NB_PHI,size_X,size_Z);
% 
% for (phi_rank=1:NB_PHI)
%     phi_rank
%     vExB_X_PR_map(:,:)=vExB_X_map_phi(phi_rank,:,:);
%     v_data=reshape(vExB_X_PR_map,NP*size_r,1);
%     vExB_X_XZ_map(:,:)=dtinterp(finesse_mesh,finesse_mesh_dtri,v_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD)';
%     vExB_X_XZ_map_phi(phi_rank,:,:)=vExB_X_XZ_map;
% end
% vExB_X_XZ_map_phi(isnan(vExB_X_XZ_map_phi))=0;
