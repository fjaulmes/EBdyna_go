z_rank=0;

gamma_values=zeros(1,100);
z_values=zeros(1,100);

for z=0.05:0.05:5;
    z_rank=z_rank+1;
    
    % z=3/4
    
    dt=0.00025;
    t_values=(dt:dt:1000);
    integrand=exp(-t_values).*t_values.^(z-1);
    
    Gamma=0;
    Gamma=sum(integrand)*dt;
    
    %for trank=1:length(t_values)
    %    integrand=exp(-t_values(trank))*t_values(trank)^(z-1);
    %    Gamma=Gamma+integrand*dt;
    %end
    
    gamma_values(z_rank)=Gamma;
    z_values(z_rank)=z;
    
end

plot(z_values,gamma_values)