function [ RMP_coils] = BS_COMPASS_RMP_load_coils(Bu_nr,Bl_nr,Au_nr,Al_nr)
%BS_AUG_RMP_load_coils Load X,Y,Z RMP coilscoordinates from data files
%   Returns struct Bu and Bl (and possibly A) for the number of coils
%   requested
%#ok<*AGROW>
%% UNIT TEST
if nargin==0
    Bu_nr=4;
    Bl_nr=4;
    Au_nr=4;
    Al_nr=4;
    
    % Add periodicy with similar coil
    check_Bu_periodic=false;    
    check_Bl_periodic=false;
    check_Au_periodic=false;
    check_Al_periodic=false;
elseif nargin==4
    % Add periodicy with similar coil
    check_Bu_periodic=false;
    check_Bl_periodic=false;
    check_Au_periodic=false;
    check_Al_periodic=false;
else
    error('Please number of Bu, Bl, Au and Al-coils')
end

%% Load files
if exist('./input/COMPASS_RMP_coil_specs','file')~=0
    % Pre-allocate
    if Bu_nr==0
        Bu=[];
    else
        Bu(Bu_nr).coil=[];
    end
    if Bl_nr==0
        Bl=[];
    else
        Bl(Bl_nr).coil=[];
    end
    if Au_nr==0
        Au=[];
    else
        Au(Au_nr).coil=[];
    end
    if Al_nr==0
        Al=[];
    else
        Al(Al_nr).coil=[];
    end
    % Load the file
    for i=1:Bu_nr
        try
            Bu(i).coil=load(['./input/COMPASS_RMP_coil_specs/Bu',num2str(i),'.asc']);  % Load matrix with R,Z,phi-coordinates
        catch
            error(strcat('Bu coil nr:',blanks(1),num2str(i),' failed to load'));
            continue
        end
        if check_Bu_periodic
            Bu(i).coil(:,3)=Bu(i).coil(:,3)-0.25*i*pi;  % Periodicy at pi/4
        end
        Bu(i).X=Bu(i).coil(:,1).*cos(Bu(i).coil(:,3));
        Bu(i).Y=Bu(i).coil(:,1).*sin(Bu(i).coil(:,3)); % No minus since of phi from the coil data runs in counterclockwise direction as seen from top
        Bu(i).Z=Bu(i).coil(:,2);
        
        % Find and remove double points
        vec=[Bu(i).X(:),Bu(i).Y(:),Bu(i).Z(:)];         % concentrate coordinates in matrix with 3 components
        dvec=diff(vec,1);                               % Find difference of previous point
        double_point=[false; all(dvec==0,2)];                    % Find which point is double
        % Delete (the second double point)
        Bu(i).X(double_point)=[];                   
        Bu(i).Y(double_point)=[];
        Bu(i).Z(double_point)=[];
        
        % Make sure the first and last point are the same (closed coil)
        if ~all(vec(1,:)==vec(end,:),2)
            Bu(i).X(end+1)=Bu(i).X(1);
            Bu(i).Y(end+1)=Bu(i).Y(1);
            Bu(i).Z(end+1)=Bu(i).Z(1);
        end
        
        Bu(i).phase=mod(mean(Bu(i).coil(:,3))-mean(Bu(1).coil(:,3)),2*pi);
        Bu(i).ref_coil=1;
       
    end
    for i=1:Bl_nr
        try
            Bl(i).coil=load(['./input/COMPASS_RMP_coil_specs/Bl',num2str(i),'.asc']); % Load matrix with R,Z,phi-coordinates
        catch
            error(strcat('Bl coil nr:',blanks(1),num2str(i),' failed to load'));
            continue
        end
        if check_Bl_periodic
            Bl(i).coil(:,3)=Bl(i).coil(:,3)-0.25*i*pi;  % Periodicy at pi/4
        end
        Bl(i).X=Bl(i).coil(:,1).*cos(Bl(i).coil(:,3));
        Bl(i).Y=Bl(i).coil(:,1).*sin(Bl(i).coil(:,3)); % No minus since of phi from the coil data runs in counterclockwise direction as seen from top
        Bl(i).Z=Bl(i).coil(:,2);
        
         % Find and remove double points
        vec=[Bl(i).X(:),Bl(i).Y(:),Bl(i).Z(:)];         % concentrate coordinates in matrix with 3 components
        dvec=diff(vec,1);                               % Find difference of previous point
        double_point=[false; all(dvec==0,2)];                    % Find which point is double
        % Delete the second duplicate point
        Bl(i).X(double_point)=[];                   
        Bl(i).Y(double_point)=[];
        Bl(i).Z(double_point)=[];
        
        % Make sure the first and last point are the same (closed coil)
        if ~all(vec(1,:)==vec(end,:),2)
            Bl(i).X(end+1)=Bl(i).X(1);
            Bl(i).Y(end+1)=Bl(i).Y(1);
            Bl(i).Z(end+1)=Bl(i).Z(1);
        end
        
        Bl(i).phase=mod(mean(Bl(i).coil(:,3))-mean(Bl(1).coil(:,3)),2*pi);
        Bl(i).ref_coil=1;
    end
    for i=1:Au_nr
        try
            Au(i).coil=load(['./input/COMPASS_RMP_coil_specs/Au',num2str(i),'.asc']); % Load matrix with R,Z,phi-coordinates
        catch
            error(strcat('Au coil nr:',blanks(1),num2str(i),' failed to load'));
            continue
        end
        if check_Au_periodic
            Au(i).coil(:,3)=Au(i).coil(:,3)-0.25*i*pi;  % Periodicy at pi/4
        end
        Au(i).X=Au(i).coil(:,1).*cos(Au(i).coil(:,3));
        Au(i).Y=Au(i).coil(:,1).*sin(Au(i).coil(:,3)); % No minus since of phi from the coil data runs in counterclockwise direction as seen from top
        Au(i).Z=Au(i).coil(:,2);
        
         % Find and remove double points
        vec=[Au(i).X(:),Au(i).Y(:),Au(i).Z(:)];         % concentrate coordinates in matrix with 3 components
        dvec=diff(vec,1);                               % Find difference of previous point
        double_point=[false; all(dvec==0,2)];                    % Find which point is double
        % Delete (the second double point)
        Au(i).X(double_point)=[];                   
        Au(i).Y(double_point)=[];
        Au(i).Z(double_point)=[];
        
        % Make sure the first and last point are the same (closed coil)
        if ~all(vec(1,:)==vec(end,:),2)
            Au(i).X(end+1)=Au(i).X(1);
            Au(i).Y(end+1)=Au(i).Y(1);
            Au(i).Z(end+1)=Au(i).Z(1);
        end
        
        Au(i).phase=mod(mean(Au(i).coil(:,3))-mean(Au(1).coil(:,3)),2*pi);
        Au(i).ref_coil=1;
                
        switch i
            case {1,2,5,6}
                Au(i).phase=mod(mean(Au(i).coil(:,3))-mean(Au(1).coil(:,3)),2*pi);
                Au(i).ref_coil=1;
            case {3,4,7,8}
                Au(i).phase=mod(mean(Au(i).coil(:,3))-mean(Au(3).coil(:,3)),2*pi);
                Au(i).ref_coil=3;
        end
    end
    for i=1:Al_nr
        try
            Al(i).coil=load(['./input/COMPASS_RMP_coil_specs/Al',num2str(i),'.asc']); % Load matrix with R,Z,phi-coordinates
        catch
            error(strcat('Al coil nr:',blanks(1),num2str(i),' failed to load'));
            continue
        end
        if check_Au_periodic
            Al(i).coil(:,3)=Al(i).coil(:,3)-0.25*i*pi;  % Periodicy at pi/4
        end
        Al(i).X=Al(i).coil(:,1).*cos(Al(i).coil(:,3));
        Al(i).Y=Al(i).coil(:,1).*sin(Al(i).coil(:,3)); % No minus since of phi from the coil data runs in counterclockwise direction as seen from top
        Al(i).Z=Al(i).coil(:,2);
        
         % Find and remove double points
        vec=[Al(i).X(:),Al(i).Y(:),Al(i).Z(:)];         % concentrate coordinates in matrix with 3 components
        dvec=diff(vec,1);                               % Find difference of previous point
        double_point=[false; all(dvec==0,2)];                    % Find which point is double
        % Delete (the second double point)
        Al(i).X(double_point)=[];                   
        Al(i).Y(double_point)=[];
        Al(i).Z(double_point)=[];
        
        % Make sure the first and last point are the same (closed coil)
        if ~all(vec(1,:)==vec(end,:),2)
            Al(i).X(end+1)=Al(i).X(1);
            Al(i).Y(end+1)=Al(i).Y(1);
            Al(i).Z(end+1)=Al(i).Z(1);
        end
        
        Al(i).phase=mod(mean(Al(i).coil(:,3))-mean(Al(1).coil(:,3)),2*pi);
        Al(i).ref_coil=1;
                
        switch i
            case {1,2,5,6}
                Al(i).phase=mod(mean(Al(i).coil(:,3))-mean(Al(1).coil(:,3)),2*pi);
                Al(i).ref_coil=1;
            case {3,4,7,8}
                Al(i).phase=mod(mean(Al(i).coil(:,3))-mean(Al(3).coil(:,3)),2*pi);
                Al(i).ref_coil=3;
        end
    end
else
    error('Maps (.asc) with COMPASS RMP coils not found')
end

%% Error if no file has been loaded
if isempty(Bu) && isempty(Bl) && isempty(Au) && isempty(Al)
    error('No data has been loaded. Did you require RMP-fields?')
end

%% UNIT TEST - plotting
if nargout==0
    delete(findall(0,'type','figure','tag','BS_single UNIT TEST'))
    hf=figure('tag','BS_single UNIT TEST');
    ha=axes('parent',hf);
    hold(ha,'on'), grid(ha,'on'), view(ha,3)
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$y$ [m]','interpreter','latex')
    zlabel(ha,'$z$ [m]','interpreter','latex')
    title(ha,'UNIT TEST loading AUG RMP coils','interpreter','latex')
    
    for i=1:Bu_nr
        h=plot3(ha,Bu(i).X,Bu(i).Y,Bu(i).Z,'r','displayname',['Bu coil ',num2str(i)]);
        if i~=1
            hasbehavior(h,'legend',false)
        end
       text(mean(Bu(i).X),mean(Bu(i).Y),mean(Bu(i).Z),[num2str(Bu(i).phase/pi,'%1.1f'),'$\pi$'],'color',get(h,'color'))
    end
    for i=1:Bl_nr        
        h=plot3(ha,Bl(i).X,Bl(i).Y,Bl(i).Z,'b','displayname',['Bl coil ',num2str(i)]);
        if i~=1
            hasbehavior(h,'legend',false)
        end
       text(mean(Bl(i).X),mean(Bl(i).Y),mean(Bl(i).Z),[num2str(Bl(i).phase/pi,'%1.1f'),'$\pi$'],'color',get(h,'color'))
    end
    for i=1:Au_nr        
        h=plot3(ha,Au(i).X,Au(i).Y,Au(i).Z,'k','displayname',['Au coil ',num2str(i)]);
        if i~=1 && i~=3
            hasbehavior(h,'legend',false)
        end
        if Au(i).ref_coil==3
            set(h,'color','g')
        end
        text(mean(Au(i).X),mean(Au(i).Y),mean(Au(i).Z),[num2str(Au(i).phase/pi,'%1.1f'),'$\pi$'],'color',get(h,'color'))
    end    
    for i=1:Al_nr        
        h=plot3(ha,Al(i).X,Al(i).Y,Al(i).Z,'k','displayname',['Al coil ',num2str(i)]);
        if i~=1 && i~=3
            hasbehavior(h,'legend',false)
        end
        if Al(i).ref_coil==3
            set(h,'color','g')
        end
        text(mean(Al(i).X),mean(Al(i).Y),mean(Al(i).Z),[num2str(Al(i).phase/pi,'%1.1f'),'$\pi$'],'color',get(h,'color'))
    end
    axis(ha,'equal')
    legend(ha,'show')
end
RMP_coils.Bu=Bu;
RMP_coils.Bl=Bl;
RMP_coils.Au=Au;
RMP_coils.Al=Al;
RMP_coils.coil_name={'Bu','Bl','Au','Al'};

end