dq_final_dr=gradient(q_final,radial_r_value_flux);

%%
%Specifying Parameters
nx=Nradial-1;               %Number of steps in space(x)
nt=20;               %Number of time steps 
dt=0.01;              %Width of each time step
dx=2/(nx-1);         %Width of space step
x=0:dx:2;            %Range of x (0,2) and specifying the grid points
u=zeros(nx,1);       %Preallocating u
un=zeros(nx,1);      %Preallocating un
vis=0.01;            %Diffusion coefficient/viscosity
beta=vis*dt/(dx*dx); %Stability criterion (0<=beta<=0.5, for explicit)
UL=1;                %Left Dirichlet B.C
UR=q_final(nx);                %Right Dirichlet B.C
UnL=dq_final_dr(1);               %Left Neumann B.C (du/dn=UnL) 
UnR=dq_final_dr(nx);               %Right Neumann B.C (du/dn=UnR) 

%%
%Initial Conditions: A square wave
for i=1:nx

        u(i)=q_final(i);

end

%%
%B.C vector
bc=zeros(nx-2,1);
bc(1)=vis*dt*UL/dx^2; bc(nx-2)=vis*dt*UR/dx^2;  %Dirichlet B.Cs
%bc(1)=-UnL*vis*dt/dx; bc(nx-2)=UnR*vis*dt/dx;  %Neumann B.Cs
%Calculating the coefficient matrix for the implicit scheme
E=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
A=E+E'-2*speye(nx-2);        %Dirichlet B.Cs
%A(1,1)=-1; A(nx-2,nx-2)=-1; %Neumann B.Cs
D=speye(nx-2)-(vis*dt/dx^2)*A;

%%
%Calculating the velocity profile for each time step
i=2:nx-1;
for it=0:nt
    un=u;
    h=plot(x,u);       %plotting the velocity profile
    axis([0 2 1 1.5])
    title({['1-D Diffusion with \nu =',num2str(vis),' and \beta = ',num2str(beta)];['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Transport property profile (u) \rightarrow')
    drawnow; 
    refreshdata(h)
    %Uncomment as necessary
    %-------------------
    %Implicit solution
    
    U=un;U(1)=[];U(end)=[];
    U=U+bc;
    U=D\U;
    u=[UL;U;UR];                      %Dirichlet
    %u=[U(1)-UnL*dx;U;U(end)+UnR*dx]; %Neumann
    %}
    %-------------------
    %Explicit method with F.D in time and C.D in space
    %{
    u(i)=un(i)+(vis*dt*(un(i+1)-2*un(i)+un(i-1))/(dx*dx));
    %}
end

for i=1:nx

        q_final_profile(i)=u(i);

end

q_final_profile(Nradial)=q_final(Nradial);

q_final_profile_diff=q_final_profile;


% close all;
figure(9)
hold on;
set(gca,'fontsize',22)
grid on
plot(radial_r_value_flux,q_final,'b--','linewidth',3)
plot(radial_r_value_flux,q_final_profile_diff,'r','linewidth',3)
xlabel('r (m)')
ylabel('q')
xlim([0 1.1])
ylim([0.8 3.4])
legend('after crash','after current diffusion')

FILENAME=strcat('../data_tokamak/q_profile.mat')
save (FILENAME,'-append','q_final_profile_diff');
