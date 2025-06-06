%% code to calculate Suzuki base cross section

function sigma=sigma_Suzuki(E_i,ne_i,Te_i,Zeff_i,A,B)
   
%sigma cm^2
%E keV/amu
%ne m^(-3)
%Te keV
%Zeff

%% test parameters(have to be commented)

% E_i=50;
% 
% ne_i=ones(50,1)*linspace(1e18,1e21,1000);
% %Te_i=ones(50,1)*linspace(0,10,1000);
% Te_i=ones(size(ne_i))*1;
% Zeff_i=ones(size(ne_i))*1;
% A=[-5.29e1 -1.36 7.19e-2 1.37e-2 4.54e-1 4.03e-1 -2.2e-1 6.66e-2 -6.77e-2 -1.48e-3];
% B(1,1,1)=-1;
% B(1,1,2)=-2.55e-2;
% B(1,2,1)=-1.25e-1;
% B(1,2,2)=-1.42e-2;
% B(2,1,1)=3.88e-1;
% B(2,1,2)=2.06e-2;
% B(2,2,1)=2.97e-2;
% B(2,2,2)=3.26e-3;
% B(3,1,1)=-2.46e-2;
% B(3,1,2)=-1.31e-3;
% B(3,2,1)=-1.48e-3;
% B(3,2,2)=-1.8e-4;


%% Decompose ne, Te, Zeff to 1d array
ne=ne_i(:);
ne(ne<1e18 & ne~=0)=1e18;
Te=Te_i(:);
Te(Te<1)=1;
Zeff=Zeff_i(:);
Zeff(Zeff<1)=1;
E=ones(size(ne))*E_i;

% 
N=ne/1e19;
lnN=log(N);
e=log(E);
U=log(Te);

for i=1:3
    for j=1:2
        for k=1:2
        Em(i,j,k,:)=e.^(i-1);
        lnNm(i,j,k,:)=lnN.^(j-1);
        Um(i,j,k,:)=U.^(k-1);
        end
    end
end
% %% 
sigH=A(1).*1e-16./E.*(1+A(2).*e+A(3).*e.^2)...
    .*(1+(1-exp(-A(4).*N)).^(A(5)).*(A(6)+A(7).*e+A(8).*e.^2))...
    .*(1+A(9).*U+A(10).*U.^2);
% 
% %sigH=1e-16/E*exp(sum(A.*Em.*lnNm.*U,'all'));
Sz=squeeze(sum(B.*Em.*lnNm.*Um,1));%summ should be not for all
Sz=squeeze(sum(Sz));
Sz=squeeze(sum(Sz));
Sz=Sz';

sigma=sigH.*(1+(Zeff-1).*Sz);
% 
% %% shape arrays back to ne,Te shape
sigma=reshape(sigma',size(ne_i));
%sigma=1;
end



