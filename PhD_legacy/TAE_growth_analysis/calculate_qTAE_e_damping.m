


qTAE=(m+0.5)/n
kpll=1/(2*R0*qTAE);

if qTAE<0.5*max(q_values)
    
    
    % Ekin_bins_lim1=0.5*(me/eV)*(vA^2+vperp_bins_lim.^2);
    % Ekin_bins_lim2=0.5*(me/eV)*((vA/3)^2+vperp_bins_lim.^2);
    % Ekin_bins1=0.5*(me/eV)*((vA/3)^2+vperp_bins_lim.^2);
    % Ekin_bins2=0.5*(me/eV)*((vA/3)^2+vperp_bins_lim.^2);
%     Ekin_bins_lim_max=max(0.5*(me/eV)*(vperp_bins_lim.^2))
%     Ekin_bins_lim=0.5*(me/eV)*(vperp_bins_lim.^2);
%     Ekin_bins=0.5*(me/eV)*(vperp_bins.^2);
    
    psi_pos_TAE=interp1(q_values,(1:length(psi_pos_bins)),qTAE)
    Ne_value=interp1(1:Nradial,n_bulk_profile,psi_pos_TAE*PSI_BIN_SIZE);
    B_value=interp1(1:Nradial,Btot_radial_profile,psi_pos_TAE*PSI_BIN_SIZE);
    vA=B_value/sqrt(mu0*Ne_value*mDT);
    wTAE=kpll*vA
    DV=0.12*vA;
    
    
%     Fpsi_pop_inf=(alphas_psi>=(psi_pos_TAE-1.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE-0.5)*PSI_BIN_SIZE);
    Fpsi_pop=(alphas_psi>=(psi_pos_TAE-0.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE+0.5)*PSI_BIN_SIZE);
    volume_slice=interp1(1:Nradial,volume_flux,(psi_pos_TAE+0.5)*PSI_BIN_SIZE)-interp1(1:Nradial,volume_flux,(psi_pos_TAE-0.5)*PSI_BIN_SIZE)
    P_slice=interp1(1:Nradial,P_profile,psi_pos_TAE*PSI_BIN_SIZE);
    n_slice=interp1(1:Nradial,n_profile,psi_pos_TAE*PSI_BIN_SIZE);
    RES_POP_SLICE=find(Fpsi_pop);
%     integrand_speed=histc(alphas_speed(RES_POP_SLICE),speed_bins_lim);
%     F_integ=0;
%     for i=1:length(speed_bins)
%         F_integ=F_integ+4*pi*(integrand_speed(i)*speed_bins(i)^2)*SPEED_BIN_SIZE;
%     end
    Npart_slice=length(RES_POP_SLICE);
%     F_integ=F_integ/Npart_slice;
%     if F_integ~=0
%         Fnorm_rescale=n_slice/F_integ;
%     else
%         Fnorm_rescale=0
%     end
    density_part_ratio=n_slice*volume_slice/Npart_slice;
%     speed_slice_sq=sum(4*pi*alphas_speed(RES_POP_SLICE).^3)/3/Npart_slice;
%     Fnorm_rescale=n_slice*volume_slice/speed_slice_sq;
    disp('Number of alphas in radial slice=')
    Npart_slice
    %%
    F_dist_vpll=zeros(length(vpll_bins),1);
    F_rescale_vpll=zeros(length(vpll_bins),1);
    integrand_speed_2D=hist2d([alphas_vpll(RES_POP_SLICE) alphas_vperp(RES_POP_SLICE)],vpll_bins_lim,speed_bins_lim);
    F_integ=0;
    for j=1:length(vpll_bins)
        F_integ_pll=0;
        for i=1:length(speed_bins)
            F_integ_pll=F_integ_pll+(integrand_speed_2D(j,i)*2*pi*speed_bins(i))*SPEED_BIN_SIZE^2;
        end
        F_integ=F_integ+F_integ_pll;
        F_dist_vpll(j)=F_integ_pll/Npart_slice;
    end
    F_dist_vpll=F_dist_vpll/sum(F_dist_vpll);
    F_rescale_vpll=n_slice.*F_dist_vpll/(VBIN_SIZE);
    
    F_integ=F_integ/Npart_slice;
    if F_integ~=0
        Fnorm_rescale_2D=n_slice/F_integ;
    else
        Fnorm_rescale_2D=0;
    end
    

%     Fpsi_pop_sup=(alphas_psi>=(psi_pos_TAE+0.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE+1.5)*PSI_BIN_SIZE);
    vpll_pop1=(alphas_vpll>=vA-DV).*(alphas_vpll<=vA+DV);
    vpll_pop2=(alphas_vpll>=vA/3-DV).*(alphas_vpll<=vA/3+DV);



    
    Fvperp1=zeros(length(vperp_bins),1);
    dFvperp1=zeros(length(vperp_bins),1);
    
    RES_POP=find(Fpsi_pop.*vpll_pop1);
    vpll_pop_size1=length(RES_POP);
    if (length(RES_POP)~= 0) && (Npart_slice~=0) && (vpll_pop_size1>3)
        % F is here per unit speed
        Frescale1=interp1(vpll_bins,F_rescale_vpll,vA)/vpll_pop_size1;
        Fvperp1=Frescale1*histc(alphas_vperp(RES_POP),vperp_bins_lim);
        dFvperp1=gradient(Fvperp1(1:end-1),vperp_bins);
    else
        Frescale1=0;
    end
    
%     FEkin1=histc(alphas_Ekin(RES_POP),Ekin_bins_lim);
%     dFEkin1=gradient(FEkin1(1:end-1),Ekin_bins);
    Fvperp2=zeros(length(vperp_bins),1);
    dFvperp2=zeros(length(vperp_bins),1);
    
    RES_POP=find(Fpsi_pop.*vpll_pop2);
    vpll_pop_size2=length(RES_POP);
    if (length(RES_POP)~= 0) && (Npart_slice~=0) && (vpll_pop_size2>3)
        % F is here per unit speed
        Frescale2=interp1(vpll_bins,F_rescale_vpll,vA/3)/vpll_pop_size2;
        Fvperp2=Frescale2*histc(alphas_vperp(RES_POP),vperp_bins_lim);
        dFvperp2=gradient(Fvperp2(1:end-1),vperp_bins);
    else
        Frescale2=0;
    end
    
%     FEkin2=histc(alphas_Ekin(RES_POP),Ekin_bins_lim); 
%     dFEkin2=gradient(FEkin2(1:end-1),Ekin_bins);
    
    Fpsi1=zeros(length(vperp_bins),length(psi_bins));
    dFpsi2=zeros(length(vperp_bins),length(psi_bins));
    Fpsi2=zeros(length(vperp_bins),length(psi_bins));
    dFpsi2=zeros(length(vperp_bins),length(psi_bins));
    

%%
    
    for (vperp_pos=VPERP0_BIN_POS:length(vperp_bins))
        
        vperp_pop=(alphas_vperp>=(vperp_pos-0.5)*VBIN_SIZE).*(alphas_vperp<=(vperp_pos+0.5)*VBIN_SIZE);
        
        %     figure(1)
        %     grid on;hold on;
        RES_POP=find(vpll_pop1.*vperp_pop);
        if (length(RES_POP)~= 0) && (Npart_slice~=0)
            Frescale1_psi=Frescale1;
            Fp1=Frescale1_psi*histc(alphas_psi(RES_POP),psi_bins_lim);
            Fpsi1(vperp_pos,:)=Fp1(1:end-1);
            dFpsi1(vperp_pos,:)=gradient(Fp1(1:end-1),psi_bins);
        else
            Fpsi1(vperp_pos,:)=zeros(length(psi_bins),1);
            dFpsi1(vperp_pos,:)=zeros(length(psi_bins),1);
        end
        RES_POP=find(vpll_pop2.*vperp_pop);
        if (length(RES_POP)~= 0) && (Npart_slice~=0)
            Frescale2_psi=Frescale2;
            Fp2=Frescale2_psi*histc(alphas_psi(RES_POP),psi_bins_lim);
            Fpsi2(vperp_pos,:)=Fp2(1:end-1);
            dFpsi2(vperp_pos,:)=gradient(Fp2(1:end-1),psi_bins);
        else
            Fpsi2(vperp_pos,:)=zeros(length(psi_bins),1);
            dFpsi2(vperp_pos,:)=zeros(length(psi_bins),1);
        end
        %     plot(q_values,Fpsi1(1:end-1));plot(q_values,Fpsi2(1:end-1),'r')
        % plot(q_values,dFpsi1);hold on;plot(q_values,dFpsi2,'r')
        
    end
    
    %%
    integ1=0;
    integ1_Ekin=0;
    integ2=0;
    integ2_Ekin=0;
    
    
    for (vperp_pos=VPERP0_BIN_POS:VPERP0_MAX_BIN_POS)
        %     vpll=vA;
        %     vperp=vperp_pos*VBIN_SIZE;
        %     integ1_Ekin=integ1_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin1(vperp_pos))*VBIN_SIZE;
        %     integ1=integ1+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin1(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi1(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
        %     vpll=vA/3;
        %     integ2=integ2+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin2(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi2(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
        %     integ2_Ekin=integ2_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin2(vperp_pos))*VBIN_SIZE;
        vperp=speed_bins(vperp_pos);
        vpll=vA;
        integ1_Ekin=integ1_Ekin+((wTAE)*(eV/me)*(vpll^2+0.5*vperp^2)*interp1(vperp_bins,dFvperp1,vperp)*SPEED_BIN_SIZE);
        integ1=integ1+((wTAE)*(eV/me)*(vpll^2+0.5*vperp^2)*interp1(vperp_bins,dFvperp1,vperp)-...
            (n)*(vpll^2+0.5*vperp^2)*vperp*interp2(vperp_bins,1:length(psi_bins),dFpsi1',vperp,psi_pos_TAE))*SPEED_BIN_SIZE;
        vpll=vA/3;
        integ2=integ2+((wTAE)*(eV/me)*(vpll^2+0.5*vperp^2)*interp1(vperp_bins,dFvperp2,vperp)-...
            (n)*(vpll^2+0.5*vperp^2)*vperp*interp2(vperp_bins,1:length(psi_bins),dFpsi2',vperp,psi_pos_TAE))*SPEED_BIN_SIZE;
        integ2_Ekin=integ2_Ekin+((wTAE)*(eV/me)*(vpll^2+0.5*vperp^2)*interp1(vperp_bins,dFvperp2,vperp))*SPEED_BIN_SIZE;
    end
    
    gamma1=(1/eV)*(1/Bavg^2)*kpll*vA*2*(pi^2)*mu0*(me^2)*R0*(qTAE^3)*integ1
    gamma2=(1/eV)*(1/Bavg^2)*kpll*vA*2*(pi^2)*mu0*(me^2)*R0*(qTAE^3)*integ2
    gamma_TAE=gamma1+gamma2
    relative_gamma=gamma_TAE/wTAE
    
    RES_POP=find(Fpsi_pop);
    disp('safety factor position')
    psi_rank=interp1(1:length(psi_pos_bins),psi_pos_bins,psi_pos_TAE)
    
    
end
