%% function to generate Ionization probobility

function [P, sigma, dP,I]=deposition_pbty(dl,E,ne,Te,Zeff,Ab,Ap,IMP,BT)
ne(isnan(ne))=0;
Te(isnan(Te))=0;
Zeff(isnan(Zeff))=0;

%% write array's A and B depending on Ap(plasma spcs atomic number) and impurity
if (E<1e2)
    if (Ap==1)
        A=[-5.29e1 -1.36 7.19e-2 1.37e-2 4.54e-1 4.03e-1 -2.2e-1 6.66e-2 -6.77e-2 -1.48e-3];
    elseif (Ap==2)
        if (BT==5)
            A=[-6.79e1 -1.22 8.14e-2 1.39e-2 4.54e-1 4.65e-1 -2.73e-1 7.51e-2 -6.3e-2 -5.08e-4];
        elseif (BT==1)
            A=[-6.98e1 -1.21 8.33e-2 1.35e-2 4.45e-1 4.89e-1 -2.90e-1 7.86e-2 -6.3e-2 -5.12e-4];  
        else 
            error('B must be aproximated by 1 or 5 [1 5]')
        end
    else
    error('wrong Ab value')
    end
elseif (E>=1e2)
    if (Ap==1)
        A=[1.27e1 1.25 4.52e-1 1.05e-2 5.47e-1 -1.02e-1 3.60e-1 -2.98e-2 -9.59e-2 4.21e-3];
    elseif (Ap==2)
        if (BT==5)
            A=[1.41e1 1.11 4.08e-1 1.05e-2 5.47e-1 -4.03e-2 3.45e-1 -2.88e-2 -9.71e-2 4.74e-3];
        elseif (BT==1)
            A=[-2.41e1 -1.30 -1.54e-1 8.02e-3 4.81e-1 -1.49e-1 3.92e-1 -2.99e-2 -9.76e-2 4.79e-3];  
        else 
            error('B must be aproximated by 1 or 5 [1 5]')
        end
    else
    error('wrong Ab value')
    end
end

if (IMP==1)%for carbon and deuterium plasmas
B(1,1,1)=-1;
B(1,1,2)=-2.55e-2;
B(1,2,1)=-1.25e-1;
B(1,2,2)=-1.42e-2;
B(2,1,1)=3.88e-1;
B(2,1,2)=2.06e-2;
B(2,2,1)=2.97e-2;
B(2,2,2)=3.26e-3;
B(3,1,1)=-2.46e-2;
B(3,1,2)=-1.31e-3;
B(3,2,1)=-1.48e-3;
B(3,2,2)=-1.8e-4;
end
%%
sigma=sigma_Suzuki(E/Ab,ne,Te,Zeff,A,B);
sigma(isnan(sigma))=0;

I=ones(size(sigma));
dI=zeros(size(sigma));

%for RK4 need to interpolate sigmas and ne to 2 time smaller grid


for i=2:(size(sigma,1)-1)
    %Euler
    %dI(i-1,:)=ne(i,:).*sigma(i,:).*I(i-1,:).*dl.*1e-4;
    %I(i,:)=I(i-1,:)-dI(i-1,:);
    %RK4
    k1=ne(i,:).*sigma(i,:).*I(i-1,:).*1e-4;
    k2=0.5*ne(i,:).*(sigma(i,:)+sigma(i+1,:)).*(I(i-1,:)+dl./2.*k1).*1e-4;
    k3=0.5*ne(i,:).*(sigma(i,:)+sigma(i+1,:)).*(I(i-1,:)+dl./2.*k2).*1e-4;
    k4=ne(i,:).*sigma(i,:).*(I(i-1,:)+dl.*k3).*1e-4;    
    Id=dl./6.*(k1+2.*k2+2.*k3+k4);
    Id(Id<0)=0;
    dI(i-1,:)=Id;
    I(i,:)=I(i-1,:)-dI(i-1,:);
end

I(I<0)=0;

I(end,:)=I(end-1,:);
dP=dI;
P=1-I;

end