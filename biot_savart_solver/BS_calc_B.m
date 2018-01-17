function [BX,BY,BZ,AX,AY,AZ] = BS_calc_B(BS,X,Y,Z,method)
%%BS_calc_B Use Biot-Savart to calculate the magnetic field. Code
%%improved from BS Toolbox 20150507 from Matlab file-exchange.
%  USE:
%    [BS,BX,BY,BZ] = BS_get_B(BS,X,Y,Z)
%
%  INPUTS:
%    BS         = BS data structure
%    X          = Field points x-coordinate vector or matrix
%    Y          = Field points y-coordinate vector or matrix
%    Z          = Field points z-coordinate vector or matrix
%
%  OUTPUTS:
%    BS         = BS data structure (no update)
%    BX         = Field points B x-component vector or matrix
%    BY         = Field points B y-component vector or matrix
%    BZ         = Field points B z-component vector or matrix

%     method=1; % Calculate B-field contribution of all points, 1 filament-element at a time
%     method=2; % Calculate B-field contribution of a coil point by point
%     method=3; % Calculate B-field contribution of all points, 1 filament-element at a time
%     method=4; % Same as method 3, but with vector potential
%     method=5; % Same as method 3, but only vector potential
%     method=6; % Make use of large matrices, but without for-loops
%     method=7; % Make use of large matrices, but without for-loops (only A)
%----------------------------------------------------
%% NARGINCHECK AND WAITBAR

if nargin~=0 && nargout~=0
    narginchk(4,5);
    if nargin==5
        method_order=method; % Try this method.
    else
        method_order=[4 3 1 2];
    end
elseif nargout==0
    % UNIT TEST
    method_order=[6 4 3 1 2];
end
waitbar_on=false;
% waitbar_on=true;

%% UNIT TEST: Wire
if nargout==0
    % UNIT TEST: 1 thin straight infinite long wire
    BS.Nfilament=1;
    BS.filament(1).Gamma=zeros(3,1e4); % number of wire segments
    BS.filament(1).Gamma(1,:)=linspace(0,1e3,size(BS.filament(1).Gamma,2)); % Length of wire.
    BS.filament(1).dGamma=1e9; % [m], maximal wire element size
    BS.filament(1).I=rand;   % [A], current
    
    % Coordinates for evaluation. Use close to wire to comparible with for infinite wire |B|=B_theta=mu0*I/2*pi*r
    x=max(BS.filament(1).Gamma(1,:))/2;     % [m], middle of wire
    x=[x x+0.01];     % [m], middle of wire
    y=linspace(0.1,1,50);                   % [m], vary y
    z=linspace(0,1,50);                     % [m], vary z
    [X,Y,Z]=meshgrid(x,y,z);                % Mesh for calculation function. (understandable output BX, BY, BZ)
    R=sqrt(Y.^2+Z.^2);                      % [m], distance from wire
end

%% Method selection
% If method has not been given (first iteration)
if ~exist('method','var') || isempty(method)
    for i=1:length(method_order)
        try
            switch method_order(i)
                case {1,2,3} % No vector potential
                    [BX,BY,BZ] = BS_calc_B(BS,X,Y,Z,method_order(i));
                    % Plot the outcome of a unit test
                    if nargout==0
                        unit_test(BS,BX,BY,BZ,X,Y,Z,R)
                    end
                case {4,5,6,7}
                    [BX,BY,BZ,AX,AY,AZ] = BS_calc_B(BS,X,Y,Z,method_order(i));
                    % Plot the outcome of a unit test
                    if nargout==0
                        if isempty(BX)
                            [BX,BY,BZ]=curl(X,Y,Z,AX,AY,AZ);
                        end
                        unit_test(BS,BX,BY,BZ,X,Y,Z,R,AX,AY,AZ)
                    end
            end
            return
        catch err
            if i==length(method_order)
                err.message=[   'Error in method: ',num2str(method_order(i)),' :',...
                                err.message];
                rethrow(err)
            else
                mess=[  'Error in method: ',num2str(method_order(i)),' :',...
                        err.message,...
                        'Continueing with method: ',num2str(method_order(i+1))];
                warning(mess)
            end
        end
    end
end

%% Size of evaluation points and evaluation point matrix
matrix_size=size(X);
numb_n=numel(X); % number of evaluation points (n)
s=cat(2,X(:),Y(:),Z(:));
clear X Y Z
s=permute(s,[1 3 2]);       

%% Pre-allocation (to add different coils)
B=zeros(numb_n,1,3);
if any(method==[4 5 6 7])
    A=zeros(numb_n,1,3);
end
%% Biot-Savart
% WAITBAR
delete(findall(0,'type','figure','tag','BS_calc_B'))
if waitbar_on
    old_text_interpreter=get(groot,'defaulttextinterpreter');
    set(groot,'defaulttextinterpreter','tex')
    wb=waitbar(0,['Filament: ',num2str(0),' @ ',num2str(0),'%'],'name','Biot-Savart RMP coils for each point','CreateCancelBtn','setappdata(gcbf,''canceling'',1)','tag','BS_calc_B');
    setappdata(wb,'canceling',0);
    
    switch method
        case {1,3}
            set(wb,'name','Biot-Savart RMP coils for filamant-element')
    end
end
for nF = 1:BS.Nfilament % Loop on each coil
    %% Coil coordinates
    % Load temporarily in shorter names
    Gamma = BS.filament(nF).Gamma;          % [m],  coordinates
    dGamma = BS.filament(nF).dGamma;        % [m],  maximum length filament length
    I_star = 1e-7*BS.filament(nF).I;        % [At], current factor in Biot-Savart law: I_star = mu0*I/(4*pi) = 4*pi*1e-7*I/(4*pi)
    dim_Gamma=find(size(Gamma)==3);         % [],   X,Y,Z dimension of coil coordinates
    
    % Check compatibility
    if ~numel(dGamma)==1; error('Maximum filament size must be a scalar'); end;
    if isempty(dim_Gamma); error('Gamma must have 1 dimension with 3-elements, indicating X,Y,Z-coordinates of coil'); end;
    
    % Change RAM ordening for optimize use of algorithm AND easy
    switch dim_Gamma
        case 1
            Gamma=permute(Gamma,[3 2 1]);       % nxlx3
        case 2
            Gamma=permute(Gamma,[3 1 2]);       % nxlx3
        otherwise
            error('Gamma must be a 2D-matrix with a 3-element dimension');
    end
	% Gamma=cat(2,Gamma,Gamma(:,1,:));  	% Ensure the coil is closed (NB: a possible duplicate point will be removed)
	
    dls=diff(Gamma,1,2);        % Difference vectors
    dls=sqrt(sum(dls.^2,3));    % Length difference vectors
    NP = ceil(dls/dGamma);      % length of filaments-elements / maximum size (i.e. cut in element size?)
    NP(NP<=0) =1;               % Interpolate all coordinates in a double occuring point
    clear dls
    
    % Now we will use index numbers for the final (maybe lengthened) filament elements
    n_a=[0 , cumsum(NP)];   %Index in final element-size of original filament coordinates (starting from 0)
    numb_l=n_a(end);        %Final number of coil coordinates
    
    p=zeros(1,numb_l+1,3);  % Coil coordinates    (Beginning / end of filament)
    l=zeros(1,numb_l  ,3);  % Filament coordinate (mid of filament)
    p(1,:,1)=interp1(n_a,Gamma(:,:,1),0:n_a(end),'linear');
    p(1,:,2)=interp1(n_a,Gamma(:,:,2),0:n_a(end),'linear');
    p(1,:,3)=interp1(n_a,Gamma(:,:,3),0:n_a(end),'linear');
    dl=diff(p,1,2);     % Difference vectors
    clear p
    
    % l-matrix, the coordinates of the filament
    l(1,:,1)=interp1(n_a,Gamma(:,:,1),0.5+0:n_a(end),'linear');
    l(1,:,2)=interp1(n_a,Gamma(:,:,2),0.5+0:n_a(end),'linear');
    l(1,:,3)=interp1(n_a,Gamma(:,:,3),0.5+0:n_a(end),'linear');
    
    switch method
        case 1
            %% Filament wise 1 
            waitbar_update=round(0.05*numb_l);  % 5 procent waitbar update
            if waitbar_on; set(wb,'name','Biot-Savart RMP coils for filamant-element'); end
            dl_x_r=zeros(numb_n,numb_l,3);     % Matrix to store cross product
            
            for m=1:numb_l
                % Difference vector
                r=s-repmat(l(:,m,:),[numb_n 1 1]);
                r=r.*repmat((sum(r.^2,3)).^(-1.5),[1 1 3]); % Normalize by |r|^3
                dl_x_r(:,m,:)=fast_cross(repmat(dl(:,m,:),[numb_n 1 1]),r,3); % Cross product
                
                % WAITBAR
                if waitbar_on && mod(m,waitbar_update)==0
                    waitbar(m/numb_l,wb,...
                        ['Filament: ',num2str(nF),' of ',num2str(BS.Nfilament),' @ ',num2str(m/numb_l*100,'%10.1f'),'%'])
                    if getappdata(wb,'canceling')
                        delete(wb)
                        error('Cancel button has been pressed')
                    end
                end
            end
            %% Store Magnetic field for this coil
            B=B+I_star*sum(dl_x_r,2); % Store vector of magnetic field for this position
        case 2
            %% Point-wise (but a whole coil applied on a point at once)
            waitbar_update=round(0.05*numb_n);
            if waitbar_on; set(wb,'name','Biot-Savart RMP coils for filamant-element'); end
            dl_x_r=zeros(numb_n,numb_l,3);     % Matrix to store cross product
            
            % Difference vector
            for n=1:numb_n
                r=repmat(s(n,:,:),[1,numb_l,1])-l;
                r=r.*repmat((sum(r.^2,3)).^(-1.5),[1 1 3]); % Normalize by |r|^3
                dl_x_r(n,:,:)=fast_cross(dl,r,3);                  % Cross product
                
                if waitbar_on && mod(n,waitbar_update)==0
                    waitbar(n/numb_n,wb,...
                        ['Coil number: ',num2str(nF),' of ',num2str(BS.Nfilament),' @ ',num2str(n/numb_n*100,'%10.1f'),'%'])
                    if getappdata(wb,'canceling')
                        delete(wb)
                        error('Cancel button has been pressed')
                    end
                end
            end
            %% Store Magnetic field for this coil
            B=B+I_star*sum(dl_x_r,2); % Store vector of magnetic field for this position
        case 3
            %% Filament wise 2 (without too large dl_x_dr-matrix)
            waitbar_update=round(0.05*numb_l);
            if waitbar_on;	set(wb,'name','Biot-Savart RMP coils for filamant-element'); end            
            B2=zeros(numb_n,1,3);
            for m=1:numb_l
                % Difference vector
                r=s-repmat(l(:,m,:),[numb_n 1 1]);
                r=r.*repmat((sum(r.^2,3)).^(-1.5),[1 1 3]); % Normalize by |r|^3
                B2=B2+fast_cross(repmat(dl(:,m,:),[numb_n 1 1]),r,3); % Cross product
                
                % WAITBAR
                if waitbar_on && mod(m,waitbar_update)==0
                    waitbar(m/numb_l,wb,...
                        ['Filament: ',num2str(nF),' of ',num2str(BS.Nfilament),' @ ',num2str(m/numb_l*100,'%10.1f'),'%'])
                    if getappdata(wb,'canceling')
                        delete(wb)
                        error('Cancel button has been pressed')
                    end
                end
            end
            %% Store Magnetic field for this coil
            B=B+B2*I_star;
        case 4
            %% Filament wise 2 (without too large dl_x_dr-matrix)
            waitbar_update=round(0.05*numb_l);
            if waitbar_on;	set(wb,'name','Biot-Savart RMP coils for filamant-element'); end            
            A2=zeros(numb_n,1,3);
            B2=zeros(numb_n,1,3);
            for m=1:numb_l
                % Difference vector
                r=bsxfun(@minus,s,l(:,m,:));
                r_norm=(sum(r.^2,3)).^(-0.5);       % Normalization: 1/|r|
                r=bsxfun(@times,r_norm.^3,r);       % Normalize by |r|^3
                
				A2=A2+bsxfun(@times,r_norm,dl(:,m,:)); % Simple vector potential
                B2=B2+fast_cross(dl(:,m,:),r,3);% Cross product
                
                % WAITBAR
                if waitbar_on && mod(m,waitbar_update)==0
                    waitbar(m/numb_l,wb,...
                        ['Filament: ',num2str(nF),' of ',num2str(BS.Nfilament),' @ ',num2str(m/numb_l*100,'%10.1f'),'%'])
                    if getappdata(wb,'canceling')
                        delete(wb)
                        error('Cancel button has been pressed')
                    end
                end
            end
            %% Store Magnetic field for this coil
            B=B+B2*I_star;
            A=A+A2*I_star;
        case 5
            %% Filament wise 2 (without too large dl_x_dr-matrix)
            waitbar_update=round(0.05*numb_l);
            if waitbar_on;	set(wb,'name','Biot-Savart RMP coils for filamant-element'); end            
            A2=zeros(numb_n,1,3);
            for m=1:numb_l
                % Difference vector
                r=bsxfun(@minus,s,l(:,m,:));
                r_norm=(sum(r.^2,3)).^(-0.5);       % Normalization: 1/|r|

                A2=A2+bsxfun(@times,r_norm,dl(:,m,:)); % Simple vector potential

                % WAITBAR
                if waitbar_on && mod(m,waitbar_update)==0
                    waitbar(m/numb_l,wb,...
                        ['Filament: ',num2str(nF),' of ',num2str(BS.Nfilament),' @ ',num2str(m/numb_l*100,'%10.1f'),'%'])
                    if getappdata(wb,'canceling')
                        delete(wb)
                        error('Cancel button has been pressed')
                    end
                end
            end
            %% Store Magnetic field for this coil
            A=A+A2*I_star;
        case 6
            %% ALL AT ONCE
            % Difference vector
            r=bsxfun(@minus,s,l);               % r=s-l
            r_norm=(sum(r.^2,3)).^(-0.5);       % Normalization: 1/|r|
            r=bsxfun(@times,r_norm.^3,r);       % Normalize by |r|^3
            
            A2 =	sum(bsxfun(@times,r_norm,dl),2); % Simple vector potential
            B2 =	sum(fast_cross(dl,r,3)      ,2);
            %% Store Magnetic field for this coil
            B=B+B2*I_star;
            A=A+A2*I_star;
        case 7
            %% ALL AT ONCE
            % Difference vector           
            A2 =sum(bsxfun(@times,(sum(bsxfun(@minus,s,l).^2,3)).^(-0.5),dl),2); % Simple vector potential
            %% Store Magnetic field for this coil
            A=A+A2*I_star;
    end
end

if waitbar_on
    delete(wb);
    set(groot,'defaulttextinterpreter',old_text_interpreter);
end

%% Reshape the output to be equal to the X,Y,Z-input size
switch method
    case {1,2,3}
        BX=reshape(B(:,:,1),matrix_size);
        BY=reshape(B(:,:,2),matrix_size);
        BZ=reshape(B(:,:,3),matrix_size);
    case {4,6}
        BX=reshape(B(:,:,1),matrix_size);
        BY=reshape(B(:,:,2),matrix_size);
        BZ=reshape(B(:,:,3),matrix_size);
        AX=reshape(A(:,:,1),matrix_size);
        AY=reshape(A(:,:,2),matrix_size);
        AZ=reshape(A(:,:,3),matrix_size);
    case {5,7}
        BX=[];
        BY=[];
        BZ=[];
        AX=reshape(A(:,:,1),matrix_size);
        AY=reshape(A(:,:,2),matrix_size);
        AZ=reshape(A(:,:,3),matrix_size);
    otherwise
        error('reshape output not programmed')
end



end

%% FAST CROSS PRODUCT
function cr=fast_cross(a,b,dim)
% No (inverse) permutation stuff, saves ~25% in cross pro
if isequal(size(a),size(b))
    switch dim
        case 1
            cr =cat(1,  a(2,:,:).*b(3,:,:)-a(3,:,:).*b(2,:,:),...
                        a(3,:,:).*b(1,:,:)-a(1,:,:).*b(3,:,:),...
                        a(1,:,:).*b(2,:,:)-a(2,:,:).*b(1,:,:));
        case 2
            cr =cat(2,  a(:,2,:).*b(:,3,:)-a(:,3,:).*b(:,2,:),...
                        a(:,3,:).*b(:,1,:)-a(:,1,:).*b(:,3,:),...
                        a(:,1,:).*b(:,2,:)-a(:,2,:).*b(:,1,:));
        case 3
            cr =cat(3,  a(:,:,2).*b(:,:,3)-a(:,:,3).*b(:,:,2),...
                        a(:,:,3).*b(:,:,1)-a(:,:,1).*b(:,:,3),...
                        a(:,:,1).*b(:,:,2)-a(:,:,2).*b(:,:,1));
        otherwise
            cr=cross(a,b,dim);
    end
elseif dim==3
            cr=cat(3,   bsxfun(@times,a(:,:,2),b(:,:,3))-bsxfun(@times,a(:,:,3),b(:,:,2)),...
                        bsxfun(@times,a(:,:,3),b(:,:,1))-bsxfun(@times,a(:,:,1),b(:,:,3)),...
                        bsxfun(@times,a(:,:,1),b(:,:,2))-bsxfun(@times,a(:,:,2),b(:,:,1)));
else
    cr=cross(a,b,dim);
end
end

%% UNIT TEST
function unit_test(BS,BX,BY,BZ,X,Y,Z,R,AX,AY,AZ)
% delete(findall(0,'type','figure','tag','Unit Test Bio-Savart integration'))
B_theory=@(r) 2e-7*BS.filament(1).I./r; % Theoretical profile of infinite wire

if nargin==11
    [BX2,BY2,BZ2]=curl(X,Y,Z,AX,AY,AZ);
end

%% B-field plot
delete(findall(0,'type','figure','tag',mfilename));
hf=figure('tag',mfilename);
ha=axes('parent',hf); hold(ha,'on')
xlabel(ha,'$r$ [m]','interpreter','latex'), ylabel(ha,'$|\vec{B}|$ [T]','interpreter','latex')
title(ha,'straight infinite long wire','interpreter','latex')

% Algorithm
B_algo=sqrt(BX(:).^2+BY(:).^2+BZ(:).^2);
matrix_info_direct=[R(:),B_algo(:)];       %Make a matrix with r and B in rising r
matrix_info_direct=sortrows(matrix_info_direct);
plot(ha,matrix_info_direct(:,1),matrix_info_direct(:,2),'+r','linewidth',3,'displayname','algorithm');
if nargin==11
    % Divergence
    B_algo=sqrt(BX2(:).^2+BY2(:).^2+BZ2(:).^2);
    matrix_info_divergence=[R(:),B_algo(:)];       %Make a matrix with r and B in rising r
    matrix_info_divergence=sortrows(matrix_info_divergence);
    plot(ha,matrix_info_divergence(:,1),matrix_info_divergence(:,2),'o','color',[0.85 0.325 0.098],'linewidth',3,'displayname','divergence');
end

% Evaluation points
h=scatter(ha,matrix_info_direct(:,1),B_theory(matrix_info_direct(:,1)),100,'k.','displayname','evaluation points');
% Theory
r=linspace(0,max(R(:)),100);
plot(ha,r,B_theory(r),'--b','linewidth',2,'displayname','theory');

hasbehavior(h,'legend',false)
legend(ha,'show')
%% Error plot
hf=figure('tag',mfilename);
ha=axes('parent',hf); hold(ha,'on')
xlabel(ha,'$r$ [m]','interpreter','latex'), ylabel(ha,'relative error [\%]','interpreter','latex')
title(ha,'straight infinite long wire','interpreter','latex')

ratio=matrix_info_direct(:,2)./B_theory(matrix_info_direct(:,1)); % Ratio of algo / theory on same points
h=plot(matrix_info_direct(:,1),100*abs(1-ratio),'.-k');
if nargin==11
    ratio=matrix_info_divergence(:,2)./B_theory(matrix_info_divergence(:,1)); % Ratio of algo / theory on same points
    h=plot(matrix_info_divergence(:,1),100*abs(1-ratio),'.-','color',[0.85 0.325 0.098],'displayname','divergence');
end

set(ha,'yscale','log')
hasbehavior(h,'legend',false)
%% Quiver plot
hf=figure('tag',mfilename);
ha=axes('parent',hf); hold(ha,'on')
xlabel(ha,'$x$ [m]','interpreter','latex'),
ylabel(ha,'$y$ [m]','interpreter','latex'),
zlabel(ha,'$z$ [m]','interpreter','latex'),
title(ha,'$\vec{B}$ [a.u.] for straight infinite long wire','interpreter','latex')
view(ha,[90 0])

% Wire and evaluation points
plot3(BS.filament.Gamma(1,:),BS.filament.Gamma(2,:),BS.filament.Gamma(3,:),'-k','linewidth',3,'displayname','wire')
scatter3(X(:),Y(:),Z(:),300,'r','marker','.','displayname','evaluation points')

% Result algorithm
B_Norm=1/(1.5*max(sqrt(BX(:).^2+BY(:).^2+BZ(:).^2))); %Normalize (automatically done otherwise), B in a.u.
quiver3(X,Y,Z,B_Norm*BX,B_Norm*BY,B_Norm*BZ,0,'-b','linewidth',3,'displayname','algorithm');   grid(ha,'on')
% Divergence B
if nargin==11
    quiver3(X,Y,Z,B_Norm*BX2,B_Norm*BY2,B_Norm*BZ2,0,'color',[0.85 0.325 0.098],'displayname','divergence');
end
% Theory
B_tot=B_theory(R);
Bx_theory=zeros(size(B_tot));
By_theory=-Z./R.*B_tot;
Bz_theory=Y./R.*B_tot;

h=quiver3(X,Y,Z,B_Norm*Bx_theory,B_Norm*By_theory,B_Norm*Bz_theory,0,'--g','displayname','theory vector');
plot3(X(:,:,1),Y(:,:,1),B_Norm*B_theory(R(:,:,1)),'-m','linewidth',3,'displayname','theory')

%     hasbehavior(h,'legend',false)
legend(ha,'show')
end