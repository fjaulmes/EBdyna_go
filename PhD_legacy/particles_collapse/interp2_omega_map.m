function interp_values=interp2_omega_map(scale_psi,scale_omega,FXZ_map,Dpsi,Domega,alphas_pos_psi,alphas_pos_omega,PART_LIST)


% alphas_pos_psi(PART_LIST)=min(alphas_pos_psi(PART_LIST),scale_psi(end));
alphas_pos_omega(PART_LIST)=min(alphas_pos_omega(PART_LIST),scale_omega(end));
% alphas_pos_psi(PART_LIST)=max(alphas_pos_psi(PART_LIST),scale_psi(1));
alphas_pos_omega(PART_LIST)=max(alphas_pos_omega(PART_LIST),scale_omega(1));
    

Xindex_precise=(alphas_pos_psi(PART_LIST)/Dpsi);
Zindex_precise=(alphas_pos_omega(PART_LIST)/Domega)+1;
%PHIindex=round(wrap(alphas_pos_phi,2)/DPHI);
Xindex_precise=min(Xindex_precise,scale_psi(end-1));


Xindex=floor(Xindex_precise);
Zindex=floor(Zindex_precise);
EXCEEDING_OMEGA=find(Zindex>(length(scale_omega)-1));

interp_x=Xindex_precise-Xindex;
interp_z=Zindex_precise-Zindex;
Zindex_p1=Zindex+1;
Zindex_p1(EXCEEDING_OMEGA)=2;
INDEX_LIST_1 = sub2ind(size(FXZ_map), Xindex, Zindex );
INDEX_LIST_2 = sub2ind(size(FXZ_map), Xindex, Zindex_p1 );
INDEX_LIST_3 = sub2ind(size(FXZ_map), Xindex+1, Zindex );
INDEX_LIST_4 = sub2ind(size(FXZ_map), Xindex+1, Zindex_p1 );




interp_values=(1-interp_x).*FXZ_map(INDEX_LIST_1).*(1-interp_z)+(interp_x).*FXZ_map(INDEX_LIST_3).*(1-interp_z)+(1-interp_x).*FXZ_map(INDEX_LIST_2).*(interp_z)+(interp_x).*FXZ_map(INDEX_LIST_4).*(interp_z);
