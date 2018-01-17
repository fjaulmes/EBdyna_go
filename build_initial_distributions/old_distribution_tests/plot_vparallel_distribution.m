    clear f_dist 
    clear vpll_range

    NumberX=400;
    DeltaX=1/(NumberX);
    

    
    % Main important values for description of the distribution
    w0=3.5*(10^6)*eV;
    v0=sqrt(2*w0/mHe);
    
    w0_pll=(1/3)*w0;
    w0_perp=(2/3)*w0;
    w0_pll=w0;
    w0_perp=w0;
    
    E0_pll=(1/3)*w0;
    E0_perp=(2/3)*w0;
    
    theta0=sqrt(2*w0/mHe);
    theta0_pll=sqrt(2*w0_pll/mHe);
    theta0_perp=sqrt(2*w0_perp/mHe);
    
    v0_pll=sqrt(2*E0_pll/mHe);


    vpll_range=(v0)*(-NumberX:NumberX)*DeltaX;
    vperp_range_sq=(v0^2-vpll_range.^2);
    vperp0=sqrt(v0^2-v0_pll^2);
    E_pll_range=0.5*mHe*vpll_range.^2;
    E_perp_range=0.5*mHe*vperp_range_sq;
    
    N0=1;
    

    f_dist=0.5*((4/3)*theta0^2-vpll_range.^2)/(theta0^3);
    f_dist=f_dist*(1/theta0);
   
    f_dist=f_dist/max(f_dist);
    f_dist=f_dist/mean(f_dist)/size(f_dist,2);
    
    f_dist_abs=f_dist(NumberX+1:end);

    
    %pos_vpll0=round(interp1(vpll_range(NumberX+1:end),(1:NumberX+1),v0_pll))
    
    Energy_parallel_range=0.01*(1:349);
    f_dist_Energy=interp1(1e-6*E_pll_range(NumberX+1:end)/eV,2*f_dist(NumberX+1:end)/max(f_dist),Energy_parallel_range);
    
    figure(3)
    plot(vpll_range,f_dist/max(f_dist));
    xlabel('v_{||} (m/s)');
    ylabel('fraction of total number of particles');
    grid on

    figure(4)
    plot(Energy_parallel_range,f_dist_Energy);
    xlabel('E_{||} (eV)');
    ylabel('fraction of total number of particles');
    grid on
    
    pos_Epll0=round(interp1(Energy_parallel_range,(1:349),E0_pll/eV*1e-6))


    % verification of energy distribution?
    % as one can see there is slightly too many particles
    % with hight velocity : probably the curvature of the sphere is not
    % well represented by our simple geometry
    disp('Energy ratio')
    sum(f_dist_Energy(1:117))

    sum(f_dist_Energy(118:349))
    