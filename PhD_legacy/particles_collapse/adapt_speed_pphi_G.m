%     alphas_pphi0_half=1.5*alphas_pphi0-0.5*alphas_pphi0_prev;
    

%     alphas_psi_value_corr=alphas_psi_value;

    v_phi_tilde=(alphas_pphi0+(ZHe)*alphas_psi_value_corr)./((mHe/eV)*alphas_Rpos_int);
%     alphas_vpll=v_X.*bX+v_Z.*bZ+v_phi.*bphi;
%     v_tilde=0.5*(alphas_vpll+v_phi_tilde./bphi);
%     v_X=v_tilde.*bX+(v_X-alphas_vpll.*bX);
%     v_Z=v_tilde.*bZ+(v_Z-alphas_vpll.*bZ);
%     v_phi=v_tilde.*bphi+(v_phi-alphas_vpll.*bphi);

%     v_phi=0.5*sign(v_phi).*abs(v_phi+v_phi_tilde);
    v_phi=0.5*(v_phi+v_phi_tilde);
%     v_phi=v_tilde;

