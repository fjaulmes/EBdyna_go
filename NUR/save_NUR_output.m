function NUR_out=save_NUR_output(FILENAME,Beam,Plasma,param,R_i,Z_i,phi_i,vR,vZ,vphi,Shth,rcoord)


%%

mu =    1.66053906660 * 1e-27;
mD =    2.01410177811 * mu;
eV =    1.602176634 * 1e-19;


v1=sqrt(2*Beam.E*(eV/mD));
v2=sqrt(2*Beam.E/2*(eV/mD));
v3=sqrt(2*Beam.E/3*(eV/mD));

vnorm=norm([vR vZ vphi]);
Enorm=0.5*(mD/eV)*vnorm^2;

Emain=find(Beam.E-Enorm<1e3);
E2=find(Beam.E/2-Enorm<1e3);
E3=find(Beam.E/3-Enorm<1e3);
vR(Emain)=vR(Emain)*v1./vnorm(Emain);
vR(E2)=vR(E2)*v2./vnorm(E2);
vR(E3)=vR(E3)*v3./vnorm(E3);
vZ(Emain)=vZ(Emain)*v1./vnorm(Emain);
vZ(E2)=vZ(E2)*v2./vnorm(E2);
vZ(E3)=vZ(E3)*v3./vnorm(E3);
vphi(Emain)=vphi(Emain)*v1./vnorm(Emain);
vphi(E2)=vphi(E2)*v2./vnorm(E2);
vphi(E3)=vphi(E3)*v3./vnorm(E3);

NUR_out=struct();
NUR_out.R_i=R_i;
NUR_out.Z_i=Z_i;
NUR_out.phi_i=phi_i;
NUR_out.vR=vR;
NUR_out.vZ=vZ;
NUR_out.vphi=vphi;
NUR_out.Shth=Shth;

INDEXES_NUR=randperm(length(NUR_out.R_i));

NUR_out.R_i=NUR_out.R_i(INDEXES_NUR);
NUR_out.Z_i=NUR_out.Z_i(INDEXES_NUR);
NUR_out.phi_i=NUR_out.phi_i(INDEXES_NUR);
NUR_out.vR=NUR_out.vR(INDEXES_NUR);
NUR_out.vZ=NUR_out.vZ(INDEXES_NUR);
NUR_out.vphi=NUR_out.vphi(INDEXES_NUR);


rcoord.X=rcoord.X(1:100:end,:);
rcoord.Y=rcoord.Y(1:100:end,:);
rcoord.Z=rcoord.Z(1:100:end,:);

rcoord.X=rcoord.X(:,INDEXES_NUR);
rcoord.Y=rcoord.Y(:,INDEXES_NUR);
rcoord.Z=rcoord.Z(:,INDEXES_NUR);



save(FILENAME,'NUR_out','Plasma','Beam','param','rcoord')

