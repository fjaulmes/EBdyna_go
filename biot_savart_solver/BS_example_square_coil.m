function [BX_square,BY_square,BZ_square]=BS_example_square_coil
%%BS_saddle_coils Calculate the magnetic field of AUG saddle coils

disp('Example using square coil!')

%% AUG dimensions roughly
AUG.R0=1.65;
AUG.a=0.5;

% Secondary coil of which field on-axis is known
Square(1).a=0.2;    % length side a
Square(1).b=0.2;    % length side b
Square(1).I = rand; % filament current [A]

% Open default figure to plot source points and field points
delete(findall(0,'type','figure','tag','example_3D_rectangle'))

% Figure for plotting coil
hf1=figure('tag','example_3D_rectangle'); ha1=axes('parent',hf1);
hold(ha1,'on'), grid(ha1,'on'), box(ha1,'on'), axis(ha1,'equal')
xlabel(ha1,'$x$ [m]','interpreter','latex'), ylabel(ha1,'$y$ [m]','interpreter','latex'), zlabel(ha1,'$z$ [m]','interpreter','latex')
view(ha1,3), axis(ha1,'tight')
axis(ha1,[-AUG.R0-AUG.a AUG.R0+AUG.a -AUG.R0-AUG.a AUG.R0+AUG.a -AUG.a AUG.a])

%% Struct BS
BS.Nfilament = 0; %Number of unit test 2 filaments so far

%% Coil coordinates
nr_p_side=100;
x_coil=[Square(1).a*ones(1,nr_p_side),...
    linspace(Square(1).a,-Square(1).a,nr_p_side),...
    -Square(1).a*ones(1,nr_p_side), ...
    linspace(-Square(1).a,Square(1).a,nr_p_side)];

y_coil=[linspace(-Square(1).b,Square(1).b,nr_p_side),...
    Square(1).b*ones(1,nr_p_side),...
    linspace(Square(1).b,-Square(1).b,nr_p_side),...
    -Square(1).b*ones(1,nr_p_side)];
z_coil=zeros(1,4*nr_p_side);
Gamma2=[x_coil;y_coil;z_coil]; %Matrix representation
Gamma2(:,[1,2,3]*nr_p_side)=[]; % Remove double points

dGamma = 2.9e-4; % filament max discretization step [m]
[BS] = BS_add_filament(BS,Gamma2,Square(1).I,dGamma,ha1);

%% Evaluation points
[X_eval,Y_eval,Z_eval]=meshgrid(0,0,unique([0 linspace(-2,2,1000)]));
[BX_square,BY_square,BZ_square] = BS_calc_B(BS,X_eval,Y_eval,Z_eval);

hf2=figure('tag','example_3D_rectangle'); ha2=axes('parent',hf2);
hold(ha2,'on'), grid(ha2,'on'), box(ha2,'on'), axis(ha2,'tight')
xlabel(ha2,'$z$ [m]','interpreter','latex'), ylabel(ha2,'$B_z$ [T]','interpreter','latex')
axis(ha2,'tight')
axis(ha2,[-AUG.R0-AUG.a AUG.R0+AUG.a -AUG.R0-AUG.a AUG.R0+AUG.a -AUG.a AUG.a])

hf3=figure('tag','example_3D_rectangle'); ha3=axes('parent',hf3);
hold(ha3,'on'), grid(ha3,'on'), box(ha3,'on'), axis(ha3,'tight')
xlabel(ha3,'$z$ [m]','interpreter','latex'), ylabel(ha3,'relative error [\%]','interpreter','latex')
axis(ha3,'tight')
axis(ha3,[-AUG.R0-AUG.a AUG.R0+AUG.a -AUG.R0-AUG.a AUG.R0+AUG.a -AUG.a AUG.a])

%% Plotting quivers
quiver3(ha1,X_eval,Y_eval,Z_eval,BX_square,BY_square,BZ_square,'color','b')
drawnow;
%% Plotting solution in z-direction
% Algorithm
plot(ha2,Z_eval(:),BZ_square(:),'-r','linewidth',3,'displayname','algorithm')

%% Theory
% WAITBAR
delete(findall(0,'type','figure','tag','Msgbox_UNIT TEST'))
msg=msgbox({'Calculating analytic solution square coil',' Please wait...'},'UNIT TEST','warn');
delete(findobj(msg,'string','OK'));
drawnow;

%        #2
%   -----2a-----
%   -          -
%#3 2b         2b  #1
%   -          -
%   -----2a-----
%        #4
% B_z is function of z only in center. Take care of different lengths
% a is now defined as being in (x) and b in (y)

% z-component important, so look only in x,y-plane-vectors (cross product)
% e.g. side #1: dl runs from -b to b, with distance sqrt(a^2+y^2+z^2).
% z-component cross product:
%       dl_x_r = dy*|r-r'| *sin(theta) = dy*a with |r-r'| distance
%       origin line in x,y-plane (z-component not considered). theta is
%       angle in x,y-plane (cylindrical coordinate) at point, so equal
%       to angle between vector dl and |r-r'|.


zs=linspace(min(Z_eval),max(Z_eval),1e4); % To provide smooth line
int=@(x,y,z) x.*(x.^2+y.^2+z.^2).^(-3/2);

% Function!
% Now for #1 and #3 + #2 #4  2 times each side
Bz_func=@(z) 2*Square.I*1e-7*integral(@(y) int(Square.a,y,z),-Square.b,Square.b)...
            +2*Square.I*1e-7*integral(@(y) int(Square.b,y,z),-Square.a,Square.a);

% B for smooth line
Bz=zeros(size(zs));
for i=1:length(zs)
    Bz(i)=Bz_func(zs(i));
end
% B to compare with evaluation points
Bz_theory=zeros(size(Z_eval));
for i=1:numel(Z_eval)
    Bz_theory(i)=Bz_func(Z_eval(i));
end
delete(msg)
plot(ha2,zs,Bz,'--b','linewidth',2,'displayname','theory');
scatter(ha2,0,4*1e-7*sqrt(2)*Square.I/Square.a,300,'.k','displayname','square analytic');
set(ha2,'ylim',[min(Bz_theory),max(Bz_theory)])

%% Plotting error
plot(ha3,Z_eval(:),100*abs(BZ_square(:)./Bz_theory(:)-1),'k.','displayname','error')
title(ha3,strcat('max error: ',num2str(100*max(abs(1-BZ_square(:)./Bz_theory(:)))),' [\%]'),'interpreter','latex')

legend(ha2,'show')
end