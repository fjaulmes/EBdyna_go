
GRIDFIT_SMOOTHNESS=0.05

FIESTA_data_X_gf=repmat(FIESTA.scale_R,1,length(FIESTA.scale_Z))-R0;
FIESTA_data_Z_gf=FIESTA.Psi*0;
zind=1;
for xind=1:length(FIESTA.scale_Z)
    FIESTA_data_Z_gf(zind:zind+length(FIESTA.scale_R)-1)=zeros(1,length(FIESTA.scale_R))+FIESTA.scale_Z(xind);
    zind=zind+length(FIESTA.scale_R);
end

gf_data=reshape(FIESTA.Psi,1,length(FIESTA.scale_R)*length(FIESTA.scale_Z));

% psi=0 at separatrix
gf_data=gf_data-FIESTA.psi_boundary;


CURRENT_DIRECTION=sign(FIESTA.psi_boundary-FIESTA.psi_axis)
% psi=1 at axis
if CURRENT_DIRECTION>0
    gf_data = -(gf_data/(FIESTA.psi_axis-FIESTA.psi_boundary)); 
else
    gf_data = (gf_data/(FIESTA.psi_axis-FIESTA.psi_boundary)); 
end


psi_XZ_map_gridfit=gridfit(FIESTA_data_X_gf,FIESTA_data_Z_gf,gf_data,X_scale,Z_scale,'smoothness',GRIDFIT_SMOOTHNESS);
% psi_XZ_map_gridfit(isnan(psi_XZ_map_gridfit)) = 1; 
psi_XZ_map_gridfit = psi_XZ_map_gridfit';

Z_Xpoint_OFFSET=abs(round(1.1*FIESTA.Z_xpoint/(FIESTA.scale_Z(2)-FIESTA.scale_Z(1))))

%%

disp('correcting normalized psi at X-point to be 1 after gridfit')
if CURRENT_DIRECTION>0
psi_n_XZ_map=(psi_XZ_map_gridfit-min(min(psi_XZ_map_gridfit(:,Z_Xpoint_OFFSET:end-Z_Xpoint_OFFSET))));
else
psi_n_XZ_map=-(psi_XZ_map_gridfit-max(max(psi_XZ_map_gridfit(:,Z_Xpoint_OFFSET:end-Z_Xpoint_OFFSET))));
end
psi_Xpoint=interp2(X_scale+R0,Z_scale,psi_n_XZ_map',FIESTA.R_xpoint,FIESTA.Z_xpoint)
psi_n_XZ_map=psi_n_XZ_map./(psi_Xpoint);


disp('The sign of current in FIESTA is positive : good for tokamak coordinates when current is clockwise !')
% gf_data=reshape(FIESTA.jphi,1,length(FIESTA.scale_R)*length(FIESTA.scale_Z));
% JPHI_XZ_map=gridfit(FIESTA_data_X_gf,FIESTA_data_Z_gf,gf_data,X_scale,Z_scale,'smoothness',GRIDFIT_SMOOTHNESS);
    [XF,YF] = meshgrid(FIESTA.scale_R-R0,FIESTA.scale_Z);
    [Xq,Yq] = meshgrid(X_scale,Z_scale);
    try
        JPHI_XZ_map=interp2(XF,YF,FIESTA.jphi',Xq,Yq,'makima');
    catch
        JPHI_XZ_map=interp2(XF,YF,FIESTA.jphi',Xq,Yq,'spline');
    end
JPHI_XZ_map=JPHI_XZ_map';

% % innacuracy induced by interpolation needs to be fixed
psi_XZ_map=psi_n_XZ_map*0;
if SIGN_CO_CURRENT_FIELD>0 
% 	psi_n_XZ_map=psi_XZ_map_gridfit-min(min(psi_XZ_map_gridfit));
    psi_XZ_map=psi_n_XZ_map*DPSI+psi_scale(1);
else
% 	psi_n_XZ_map=psi_XZ_map_gridfit-max(max(psi_XZ_map_gridfit))+1;
    psi_XZ_map=psi_n_XZ_map*DPSI+psi_scale(1);
    warning('check psi_XZ_map!')
end

psi_n_XZ_map_sat=min(psi_n_XZ_map,1);
radial_XZ_map=psi_n_XZ_map_sat*(Nradial-1)+1;

P_XZ_map=interp1(psi_n_scale,P_initial_profile,psi_n_XZ_map_sat);

disp('Toroidal sign convention of direct <> >0 incompatible with (R,Z,phi) tokamak coordinates')
% disp('The sign of Bphi from EFIT needs to be reversed!!!')
% this is so true for the F function !
F_XZ_map=interp1(psi_n_scale,Fdia_profile,psi_n_XZ_map_sat);
Bphi_XZ_map=F_XZ_map./Rpos_XZ_map;

derpsiX=psi_XZ_map*0;
derpsiZ=psi_XZ_map*0;

for z=3:length(Z_scale)-2
    for x=3:length(X_scale)-2
        derpsiX(x,z)=(1/12)*(-psi_XZ_map(x+2,z)+psi_XZ_map(x-2,z))+(2/3)*(psi_XZ_map(x+1,z)-psi_XZ_map(x-1,z));
        derpsiZ(x,z)=(1/12)*(-psi_XZ_map(x,z+2)+psi_XZ_map(x,z-2))+(2/3)*(psi_XZ_map(x,z+1)-psi_XZ_map(x,z-1));
        derpsiX(x,z)=derpsiX(x,z)/DX;
        derpsiZ(x,z)=derpsiZ(x,z)/DX;       
    end
end

BR_XZ_map=-derpsiZ./Rpos_XZ_map;
BZ_XZ_map=derpsiX./Rpos_XZ_map;
Bpol_XZ_map=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);


q_XZ_map=interp1(psi_n_scale,q_initial_profile,psi_n_XZ_map_sat);

[value X2_out]=max(psi_n_XZ_map_sat(mid_X:end,Z_axis_pos));
[value X1_out]=max(psi_n_XZ_map_sat(mid_X:-1:1,Z_axis_pos));

X2_out=X2_out+mid_X-1
X1_out=mid_X-X1_out+1


radial_XZ_map=interp1(psi_n_FIESTA,radial_r_values,psi_n_XZ_map_sat);

