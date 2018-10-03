
z=2/4

dt=0.00025
t_values=(dt:dt:1000);
integrand=exp(-t_values).*t_values.^(z-1);

Gamma=0;
Gamma=sum(integrand)*dt;

%for trank=1:length(t_values)
%    integrand=exp(-t_values(trank))*t_values(trank)^(z-1);
%    Gamma=Gamma+integrand*dt;
%end

gamma_value_34=Gamma


z=5/4


integrand=exp(-t_values).*t_values.^(z-1);

Gamma=0;
Gamma=sum(integrand)*dt;

gamma_value_54=Gamma

C0=2*gamma_value_54/gamma_value_34

product=gamma_value_54*gamma_value_34
product_verif=(pi/4)/sin(pi/4)