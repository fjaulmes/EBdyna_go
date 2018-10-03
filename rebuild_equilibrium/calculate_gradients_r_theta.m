
grad_r_PR_map=zeros(NP,Nradial);
grad_theta_PR_map=zeros(NP,Nradial);
grad_theta_PR_map_X=zeros(NP,Nradial);
grad_theta_PR_map_Z=zeros(NP,Nradial);

for(p=1:NP)
    for(r=2:Nradial-1)        
        dr_X=X_PR_map(p,r+1)-X_PR_map(p,r);
        dr_Z=Z_PR_map(p,r+1)-Z_PR_map(p,r);
        dr1=sqrt(dr_X^2+dr_Z^2);
        r1=sqrt(X_PR_map(p,r+1)^2+Z_PR_map(p,r+1)^2);
        
        dr_X=X_PR_map(p,r)-X_PR_map(p,r-1);
        dr_Z=Z_PR_map(p,r)-Z_PR_map(p,r-1);
        dr2=sqrt(dr_X^2+dr_Z^2);
        r2=sqrt(X_PR_map(p,r-1)^2+Z_PR_map(p,r-1)^2);
        
        grad_r_PR_map(p,r)=(r1-r2)/(dr1+dr2);
    end
end
grad_r_PR_map(:,1)=zeros(NP,1);
grad_r_PR_map(:,Nradial)=zeros(NP,1);

dl_X=X_PR_map*0;
dl_Z=X_PR_map*0;


for(p=2:NP-1)
    for(r=1:Nradial)        
        dl_X1=X_PR_map(p+1,r)-X_PR_map(p,r);
        dl_Z1=Z_PR_map(p+1,r)-Z_PR_map(p,r);
        dl1=sqrt(dl_X1^2+dl_Z1^2);
        l1=2*pi/(NP);
        
        dl_X2=X_PR_map(p,r)-X_PR_map(p-1,r);
        dl_Z2=Z_PR_map(p,r)-Z_PR_map(p-1,r);
        dl2=sqrt(dl_X2^2+dl_Z2^2);
        l2=2*pi/(NP);
        dl_X(p,r)=0.5*(dl_X1+dl_X2);
		dl_Z(p,r)=0.5*(dl_Z1+dl_Z2);
		
        grad_theta_PR_map(p,r)=(l1+l2)/(dl1+dl2);
        grad_theta_PR_map_X(p,r)=(dl_X1+dl_X2)/(dl1+dl2);
        grad_theta_PR_map_Z(p,r)=(dl_Z1+dl_Z2)/(dl1+dl2);
    end
end

for(r=1:Nradial)
    dl_X1=X_PR_map(2,r)-X_PR_map(1,r);
    dl_Z1=Z_PR_map(2,r)-Z_PR_map(1,r);
    dl1=sqrt(dl_X1^2+dl_Z1^2);
    l1=2*pi/(NP);
    
    dl_X2=X_PR_map(1,r)-X_PR_map(NP-1,r);
    dl_Z2=Z_PR_map(1,r)-Z_PR_map(NP-1,r);
    dl2=sqrt(dl_X2^2+dl_Z2^2);
    l2=2*pi/(NP);
    dl_X(1,r)=0.5*(dl_X1+dl_X2);
	dl_Z(1,r)=0.5*(dl_Z1+dl_Z2);
		
    grad_theta_PR_map(1,r)=(l1+l2)/(dl1+dl2);
    grad_theta_PR_map_X(1,r)=(dl_X1+dl_X2)/(dl1+dl2);
    grad_theta_PR_map_Z(1,r)=(dl_Z1+dl_Z2)/(dl1+dl2);
end



grad_theta_PR_map(NP,:)=grad_theta_PR_map(1,:);
grad_theta_PR_map_X(NP,:)=grad_theta_PR_map_X(1,:);
grad_theta_PR_map_Z(NP,:)=grad_theta_PR_map_Z(1,:);

grad_theta_PR_map_X=grad_theta_PR_map_X./sqrt(grad_theta_PR_map_X.^2+grad_theta_PR_map_Z.^2);
grad_theta_PR_map_Z=grad_theta_PR_map_Z./sqrt(grad_theta_PR_map_X.^2+grad_theta_PR_map_Z.^2);
grad_theta_PR_map_X=grad_theta_PR_map_X./sqrt(grad_theta_PR_map_X.^2+grad_theta_PR_map_Z.^2);
grad_theta_PR_map_Z=grad_theta_PR_map_Z./sqrt(grad_theta_PR_map_X.^2+grad_theta_PR_map_Z.^2);

sqrtg_PR_map=Rpos_PR_map./sqrt((grad_r_PR_map.^2).*grad_theta_PR_map.^2-(cos_ki_PR_map.^2).*(grad_theta_PR_map.*grad_r_PR_map).^2);
sqrtg_PR_map(:,1)=sqrtg_PR_map(:,2);
sqrtg_PR_map(:,Nradial)=sqrtg_PR_map(:,Nradial-1);


g_data=reshape((grad_theta_PR_map_X(:,:)),NP*Nradial,1);
g_polX_XZ_map=griddata(finesse_data_X,finesse_data_Z,g_data,XX,ZZ,'linear');
g_polX_XZ_map(isnan(g_polX_XZ_map)) = 0; 
g_polX_XZ_map=g_polX_XZ_map';


g_data=reshape((grad_theta_PR_map_Z(:,:)),NP*Nradial,1);
g_polZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,g_data,XX,ZZ,'linear');
g_polZ_XZ_map(isnan(g_polZ_XZ_map)) = 0; 
g_polZ_XZ_map=g_polZ_XZ_map';