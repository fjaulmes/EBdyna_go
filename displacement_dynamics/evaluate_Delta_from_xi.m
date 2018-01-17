%%
% representation of the circles
recalculate_layer_width;

DX_lin=0.0005;
Xlinscale=(0:DX_lin:MAX_DELTA_RADIUS);
acoef=zeros(length(Xlinscale),1);

deltaarray=zeros(length(Xlinscale),1);
delta_posx=zeros(length(r_core),1);
delta_acoef=zeros(length(r_core),1);
delta_theta0=zeros(1,length(r_core));
delta_XI2=zeros(1,length(r_core));
delta_XI1=zeros(1,length(r_core));
delta_YI1=zeros(1,length(r_core));

linearray=zeros(1,length(Xlinscale));
Ycore=zeros(length(Xlinscale),length(r_core));
Yrx=zeros(length(Xlinscale),length(r_core));
for t=1:length(r_core)
    Ycore(:,t)=sqrt(r_core(t)^2-Xlinscale.^2)-r_core(t);
    Yrx(:,t)=sqrt(rx_evol_lin(t)^2-Xlinscale.^2)-rx_evol_lin(t);
end
Yrx(imag(Yrx(:,:))~=0)=0;
Ycore(imag(Ycore(:,:))~=0)=0;
% Yrx(~isreal(Yrx))=0;
% Ycore(~isreal(Ycore))=0;
for t=2:TRANSITION_FRAME
    acoef(2:end-1)=(Ycore(2:end-1,t)+r_core(t))./Xlinscale(2:end-1)';
    XI1=zeros(length(Xlinscale),1);
    YI1=zeros(length(Xlinscale),1);
    XI2=zeros(length(Xlinscale),1);
    YI2=zeros(length(Xlinscale),1);
    length_core=floor(r_core(t)/DX_lin);
    for x=2:length_core
        if acoef(x)~=0
            linearray=acoef(x)*Xlinscale-rx_evol_lin(t);
            linearray=linearray';
            XI1(x)=interp1((Ycore(1:length_core,t)-linearray(1:length_core)),Xlinscale(1:length_core),0,'cubic');
            YI1(x)=interp1(Xlinscale,Ycore(:,t),XI1(x));
            XI2(x)=interp1((Yrx(1:length_core,t)-linearray(1:length_core)),Xlinscale(1:length_core),0,'cubic');
            YI2(x)=interp1(Xlinscale,Yrx(:,t),XI2(x));
        end
    end
    deltaarray=sqrt((YI2-YI1).^2+(XI2-XI1).^2);
    [valmin deltaindex]=min(abs(delta_evol(t)-deltaarray(:)));
    delta_posx(t)=deltaindex;
    delta_acoef(t)=acoef(deltaindex);
    delta_theta0(t)=pi/2-atan(delta_acoef(t));
    delta_XI2(t)=XI2(deltaindex);
    delta_XI1(t)=XI1(deltaindex);
    delta_YI1(t)=r_core(t)+YI1(deltaindex);
end
%fix first point
delta_theta0(1)=delta_theta0(2);



thetac_evol=pi/2-atan(delta_YI1(1:length(vA_out))./delta_XI1(1:length(vA_out)));
thetac_evol(1)=thetac_evol(2);

Delta_SP=2*(thetac_evol(1:length(vA_out))).*(0.5*(r_core(1:length(vA_out))+rx_evol_lin(1:length(vA_out))));
% Delta_SP=2*(delta_theta0(1:length(vA_out))).*rx_recalc_evol(1:length(vA_out));
tau_star=2*Delta_SP./vA_out;

% recalculate_layer_width;

%this should be a function of xi

%g_xi=inv_tauA.^2+mu_delta*inv_tauA;
g_xi=ksi_dot_recalc(1:TRANSITION_FRAME);
% g_xi=ksi_dot_recalc(1:TRANSITION_FRAME)*crash_duration/crash_duration_ini;
g_xi_error=(2*vA_out.*delta_evol./Delta_SP)-ksi_dot_recalc(1:TRANSITION_FRAME);
g_xi_error_integ=g_xi_error_integ+0.05*g_xi_error;
g_xi=g_xi+0.1*g_xi_error+0.1*g_xi_error_integ;

% g_xi=(2*vA_out.*delta_evol./Delta_SP);
g_xi(isnan(g_xi))=max(g_xi);
% correct the last points
g_xi(end-1)=1.5*g_xi(end-2)-0.5*g_xi(end-2);
g_xi(end)=g_xi(end-1);
