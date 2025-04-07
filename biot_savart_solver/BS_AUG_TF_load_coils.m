function [ TF] = BS_AUG_TF_load_coils
%BS_AUG_TFR_load_coils Load X,Y,Z TFR coils coordinates from data file
%   Return TF_coil data

%% Load files
if ~exist('./input/modcoil_1coil.dat','file')
    error('Map with AUG TF coil not found')
else
    % Pre-allocate
    TF(16).coil=[];
    % Load single coil data
    M=importdata('./input/modcoil_1coil.dat',' ',10); 
    % Please note the data consists of 1 rectangular shaped D-coil. So all
    % 4 corners are given. To get a single filament wire, we model
    % this by finding the middle coordinate of the coil and the current
    % through this. First we need to distinguish the 4 sides
    
%     x_TF=M.data(:,1);
%     y_TF=M.data(:,2);
%     z_TF=M.data(:,3);
    pos_same=(all(bsxfun(@eq,M.data(:,1:3),M.data(1              ,1:3)),2));                  % Find the position with all 3 coordinates equal to the begin of the 1st side
    if sum(pos_same)~=2 ; error('Not a single equal point found to close side of coil data'); end;
    ind_same(1)=find(pos_same,1,'last');
    pos_same=(all(bsxfun(@eq,M.data(:,1:3),M.data(ind_same(end)+1,1:3)),2));    % Find the position with all 3 coordinates equal to the begin of the 2nd side
    if sum(pos_same)~=2 ; error('Not a single equal point found to close side of coil data'); end;
    ind_same(2)=find(pos_same,1,'last');
    pos_same=(all(bsxfun(@eq,M.data(:,1:3),M.data(ind_same(end)+1,1:3)),2));    % Find the position with all 3 coordinates equal to the begin of the 3rd side
    if sum(pos_same)~=2 ; error('Not a single equal point found to close side of coil data'); end;
    ind_same(3)=find(pos_same,1,'last');
    pos_same=(all(bsxfun(@eq,M.data(:,1:3),M.data(ind_same(end)+1,1:3)),2));    % Find the position with all 3 coordinates equal to the begin of the 4th side
    if sum(pos_same)~=2 ; error('Not a single equal point found to close side of coil data'); end;
    ind_same(4)=find(pos_same,1,'last');
    
    points_per_side=diff([0 ind_same]);
    if any(points_per_side~=points_per_side(1))
        error('Failed to distinguish 4 sides of TF coil data')
    end
    
    new_M=zeros([size(M.data(:,1:3),1)/4 3 4]); % Reshape M to house 4 points per position (last dimension)
    new_M(1:points_per_side,:,1)=M.data(1               :   ind_same(1),1:3);
    new_M(1:points_per_side,:,2)=M.data(ind_same(1)+1   :   ind_same(2),1:3);
    new_M(1:points_per_side,:,3)=M.data(ind_same(2)+1   :   ind_same(3),1:3);
    new_M(1:points_per_side,:,4)=M.data(ind_same(3)+1   :   ind_same(4),1:3);
    
    coord_center=mean(new_M,3);
    coord_center(:,2)=0;        % Define y=0 for basis coil, x -> R, z-> Z (so now coord_center has cylindrical coordinates)
    phi_coils=linspace(0,2*pi,17); phi_coils(end)=[];
    
    x_coils= coord_center(:,1)*cos(phi_coils);
    y_coils=-coord_center(:,1)*sin(phi_coils);  % Minus since of order / right-handed axis / phi runs in negative Y
    z_coils= coord_center(:,3)*ones(size(phi_coils));
    
    % Put the data in coil struct
    for i=1:16    
        TF(i).X=x_coils(:,i);
        TF(i).Y=y_coils(:,i);
        TF(i).Z=z_coils(:,i);
    end
end

%% UNIT TEST - plotting
if nargout==0
    delete(findall(0,'type','figure','tag',mfilename))
    hf=figure('tag',mfilename);
    ha=axes('parent',hf);
    hold(ha,'on'), grid(ha,'on'), view(ha,3)
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$y$ [m]','interpreter','latex')
    zlabel(ha,'$z$ [m]','interpreter','latex')
    title(ha,'UNIT TEST loading AUG TF coils','interpreter','latex')
    
    for i=1:16
        h=plot3(ha,TF(i).X,TF(i).Y,TF(i).Z,'k','displayname',['TF coil ',num2str(i)]);
        if i~=1
            hasbehavior(h,'legend',false)
        end
    end
    axis(ha,'equal')
    legend(ha,'show')
end

end