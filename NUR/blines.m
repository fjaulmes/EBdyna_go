%% EBDYNA beam


function [X,Z]=blines(bmwid,focl,div,nb,app)
%% test_parameters
% div.r = 0.00823688265;                        % ion source horizontal divergence (radians)
% div.z = 0.00823688265;                        % ion source horizontal divergence (radians)
% bmwid.r = 0.1;                      % ion source half-width (m)
% bmwid.z = 0.2;                      % ion source half-width (m)
% bmwid.opt = 1;                      % ion source half-width (m)
% focl.r = 4;                     % ion source horizontal focal length (m)
% focl.z = 4;                     % ion source horizontal focal length (m)
% app.l(1)=3;                 % app aperture to source distance
% app.l(2)=3.2;                 % app aperture to source distance
% app.r(1)=0.08;                  % app(2) apperture radius
% app.r(2)=0.088;                  % app(2) apperture radius
% app.z(1)=0.089;                  % app(2) apperture radius
% app.z(2)=0.08;                  % app(2) apperture radius
% app.opt(1)=0;                   % option: 0 -circ, 1 - rect
% app.opt(2)=1;                   % option: 0 -circ, 1 - rect
% nb=100;                        %nb numbear of generated beams


%% generate_ ion_ birth_ position
i=0;

if (bmwid.opt==0)
    Xb=zeros(1,2*nb);
    Zb=zeros(1,2*nb);
    while(i<2*nb)
        Xbt=(1-2*rand)*bmwid.r;
        Zbt=(1-2*rand)*bmwid.r;
        if (hypot(Xbt,Zbt)<bmwid.r)
            i=i+1;
            Xb(i)=Xbt;
            Zb(i)=Zbt;
        end
    end
elseif (bmwid.opt==1)
    Xb=(1-2*rand(1,2*nb))*bmwid.r;
    Zb=(1-2*rand(1,2*nb))*bmwid.z;       
end
    

%% find _an angle_ to_ reach _the_ focal_ point_ and_ deviation_ with_ Gaussian_ distribution
Xtg=-Xb/focl.r;
Ztg=-Zb/focl.z;

%add random normally distributed value
Xtg=Xtg+normrnd(0,div.r,size(Xtg));
Ztg=Ztg+normrnd(0,div.z,size(Ztg));

%% generate_ path_ coordinates
l=[0,app.l];

X=ones(size(l'))*Xb+(l'*ones(size(Xtg))).*(ones(size(l'))*Xtg);
Z=ones(size(l'))*Zb+(l'*ones(size(Ztg))).*(ones(size(l'))*Ztg);


%% find_ beam _bad _indexes (beams that heat apperture)
% radius
R=hypot(X,Z);
ind=ones(size(Xb));
        
for i=1:numel(app.l)
    if (app.opt(i)==0)
        ind(find(R(i+1,:)>app.r(i)))=NaN;
    elseif (app.opt(i)==1)
        ind(find(abs(X(i+1,:))>app.r(i)))=NaN;
        ind(find(abs(Z(i+1,:))>app.z(i)))=NaN;
    end
end

%% remove_ bad _l,X,Z
X(:,isnan(ind))=[];
Z(:,isnan(ind))=[];
%% reduce array to nb dimentions
if (size(X,2)<nb)
    error('more then 50% of beam deposed on the apperture, DO SOMETHING! YOU HAVE BAD INJECTOR :)')
end

X=X(:,1:nb);
Z=Z(:,1:nb);

X=X(1:2,:);
Z=Z(1:2,:);

end