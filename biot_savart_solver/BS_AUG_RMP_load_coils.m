function [ RMP_coils] = BS_AUG_RMP_load_coils(Bu_nr,Bl_nr,A_nr)
%BS_AUG_RMP_load_coils Load X,Y,Z RMP coilscoordinates from data files
%   Returns struct Bu and Bl (and possibly A) for the number of coils
%   requested
%#ok<*AGROW>
%% UNIT TEST
if nargin==0
    Bu_nr=8;
    Bl_nr=8;
    A_nr=0;
    
    % Add periodicy with similar coil
    check_Bu_periodic=false;    
    check_Bl_periodic=false;
    check_A_periodic=false;
elseif nargin==3
    % Add periodicy with similar coil
    check_Bu_periodic=false;
    check_Bl_periodic=false;
    check_A_periodic=false;
else
    error('Please number of Bu, Bl and A-coils')
end

%% Load files
if exist('./input/AUG RMP coil specs','file')~=0
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
    if A_nr==0
        A=[];
    else
        A(A_nr).coil=[];
    end
    % Load the file
    for i=1:Bu_nr
        try
            Bu(i).coil=load(['./input/AUG RMP coil specs/Bu',num2str(i),'.asc']);  % Load matrix with R,Z,phi-coordinates
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
            Bl(i).coil=load(['./input/AUG RMP coil specs/Bl',num2str(i),'.asc']); % Load matrix with R,Z,phi-coordinates
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
    for i=1:A_nr
        try
            A(i).coil=load(['./input/AUG RMP coil specs/A',num2str(i),'.asc']); % Load matrix with R,Z,phi-coordinates
        catch
            error(strcat('A coil nr:',blanks(1),num2str(i),' failed to load'));
            continue
        end
        if check_A_periodic
            A(i).coil(:,3)=A(i).coil(:,3)-0.25*i*pi;  % Periodicy at pi/4
        end
        A(i).X=A(i).coil(:,1).*cos(A(i).coil(:,3));
        A(i).Y=A(i).coil(:,1).*sin(A(i).coil(:,3)); % No minus since of phi from the coil data runs in counterclockwise direction as seen from top
        A(i).Z=A(i).coil(:,2);
        
         % Find and remove double points
        vec=[A(i).X(:),A(i).Y(:),A(i).Z(:)];         % concentrate coordinates in matrix with 3 components
        dvec=diff(vec,1);                               % Find difference of previous point
        double_point=[false; all(dvec==0,2)];                    % Find which point is double
        % Delete (the second double point)
        A(i).X(double_point)=[];                   
        A(i).Y(double_point)=[];
        A(i).Z(double_point)=[];
        
        % Make sure the first and last point are the same (closed coil)
        if ~all(vec(1,:)==vec(end,:),2)
            A(i).X(end+1)=A(i).X(1);
            A(i).Y(end+1)=A(i).Y(1);
            A(i).Z(end+1)=A(i).Z(1);
        end
        
        A(i).phase=mod(mean(A(i).coil(:,3))-mean(A(1).coil(:,3)),2*pi);
        A(i).ref_coil=1;
                
        switch i
            case {1,2,5,6}
                A(i).phase=mod(mean(A(i).coil(:,3))-mean(A(1).coil(:,3)),2*pi);
                A(i).ref_coil=1;
            case {3,4,7,8}
                A(i).phase=mod(mean(A(i).coil(:,3))-mean(A(3).coil(:,3)),2*pi);
                A(i).ref_coil=3;
        end
    end
else
    error('Map with AUG RMP coils not found')
end

%% Error if no file has been loaded
if isempty(Bu) && isempty(Bl) && isempty(A)
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
    for i=1:A_nr        
        h=plot3(ha,A(i).X,A(i).Y,A(i).Z,'k','displayname',['A coil ',num2str(i)]);
        if i~=1 && i~=3
            hasbehavior(h,'legend',false)
        end
        if A(i).ref_coil==3
            set(h,'color','g')
        end
        text(mean(A(i).X),mean(A(i).Y),mean(A(i).Z),[num2str(A(i).phase/pi,'%1.1f'),'$\pi$'],'color',get(h,'color'))
    end
    axis(ha,'equal')
    legend(ha,'show')
end
RMP_coils.Bu=Bu;
RMP_coils.Bl=Bl;
RMP_coils.A=A;
RMP_coils.coil_name={'Bu','Bl','A'};

end