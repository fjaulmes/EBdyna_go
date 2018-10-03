%% Default after loading data
if ~exist('par','var') || isempty(par)
    error('parameter information in par-struct missing')
end
delete(findall(0,'type','figure','tag',mfilename));
clear hf ha h hl

title_st={['sim \# ',par.ID],['$\Delta t=',num2str(par.dt,'%1.1e'),'$ [s], $t_{sim}=',num2str(par.time_scale(end),'%1.0e'),'$ [s]']};

%% Make some variables global
make_global={'par','maps','dim','const','time','qom'};
for i=1:length(make_global)
    if ~exist(make_global{i},'var')
        eval(['global ',make_global{i}]);
    else
        a=whos(make_global{i});
        if size(a,1)==1 && ~a.global
            eval(['temp=',make_global{i},';']);
            clear(make_global{i});
            eval(['global ',make_global{i}]);
            eval([make_global{i},' =temp;']);
        end
        clear temp
    end
end

%% Select plot routine
% 1 - Slideshow of each particles orbit and physical quantities
% 2 - Phase space plot with free to choose vertical / horizontal
% 3 - Scatter plot of trapped particles wb/wd and Delta pphi
% 4 - 3D histogram of (normalized) trapped particles wb/wd
% 5 - Figures in pphi - mu space with pphi the color (expr=dpphi_dt>0.00, color up to 0.04)
% 6 - Figure in mu - q_avg space to categorize
% 7 - Figure in mu - pphi space to categorize
% 8 - Contourplots in mu - pphi space with losses

plot_routine=8; % Define plot routine

%% Find q initial and the q on the avgerage flux surface
if any(plot_routine==[1]) && (~isfield(output,'pphi_kin') || isempty(dim))
    reset_data_analysis_environment;
end
if ~exist('q_initial','var') && any(plot_routine==[1 2])
    if ~isfield(dim,'q_initial_profile')
        dim2=load([par.paths.DATA_FOLDER,'q_profile.mat'],'q_initial_profile');        dim=combine_structs(dim,dim2);
        dim2=load('../data_tokamak/motions_map_dimensions.mat','scale_X','scale_Z','R0');        dim=combine_structs(dim,dim2);
        dim2=load('../data_tokamak/psi_profiles.mat','psi_pol_initial_profile');        dim=combine_structs(dim,dim2);
        maps2=load('../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat','psi_norm_XZsmall_map'); maps.psi_norm_XZ=maps2.psi_norm_XZsmall_map;
    end
    radial_initial=interp2(dim.scale_X,dim.scale_Z,maps(1).psi_norm_XZ',input.x_gc(:,1)-dim.R0,input.x_gc(:,2));
    q_initial=interp1(1:513,dim.q_initial_profile,radial_initial);
    if exist('prec','var')
        prec.q_on_avg_FS=interp1(dim.psi_pol_initial_profile,dim.q_initial_profile,prec.psi_avg);
    end
end

a=load('../data_tokamak/B_fields.mat','Bpol_PR_map','Btor_PR_map');
B0=mean(sqrt(a.Bpol_PR_map(:,1).^2+a.Btor_PR_map(:,1).^2)); clear a


%% Find B-field at input
if (~exist('B','var') || ~exist('b','var')) && any(plot_routine==[1])    
    %% Find B (2D)
    time=0; % Put global time to zero (to get the equilibrium in ST simulations)
    [~,B] = B_interpolation(input.x);
    Bfield_sq=dot(B,B,2);
    Bfield=sqrt(Bfield_sq); clear Bfield_sq
    
    b=bsxfun(@times,Bfield.^-1,B);
end
%% Find Eperp
if ~exist('Eperp','var') && any(plot_routine==[1])
    vpll=dot(b,input.v,2);
    Eparralel=(0.5*input.m/const.eV)*vpll.^2;
    Eperp=input.Ekin-(0.5*input.m/const.eV)*vpll.^2;
end

%% PLOTTING
switch plot_routine
    case 1
        %% Movie of orbit each particle
        % Figure
        hf=figure('tag',mfilename);
        ha(1)=subplot(3,2,1,'parent',hf);
        ha(2)=subplot(3,2,3,'parent',hf);
        ha(3)=subplot(3,2,5,'parent',hf);
        ha(4)=subplot(1,2,2,'parent',hf); axis(ha(4),'equal');
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        set(ha,'FontSize',28)
        xlabel(ha(3),'$t$ [s]','FontSize',20)
        xlabel(ha(4),'$R$ [m]','FontSize',20)
        ylabel(ha(1),'$\Delta p_\varphi/\psi_0$','FontSize',20)
        ylabel(ha(2),'$\mu$ [eV/T]','FontSize',20)
        ylabel(ha(3),'$\Delta E_{kin}/E_{kin}$','FontSize',20)
        ylabel(ha(4),'$Z$ [m]','FontSize',20)
        hold(ha(1),'on') ; hold(ha(2),'on'); hold(ha(3),'on');hold(ha(4),'on')
        
        xlim(ha(1),[0 par.time_scale(end)])
        linkaxes(ha(1:3),'x');
        
        indexes=1:input.N_job;
        c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps(1).psi_norm_XZ',[dim.NB_PSI dim.NB_PSI]);
        plot(c(1,2:end),(c(2,2:end)),'r','displayname','LCFS','parent',ha(4))
        for ind=1600:1800
            % Left upper
            h(1)=plot(ha(1),par.time_scale,bsxfun(@minus,output.pphi_kin(ind,:),output.pphi_kin(ind,1))','b','displayname','$mRv_\varphi + qRA_\varphi$');
            if par.CALCULATE_PPHI_3D
                h(end+1)=plot(ha(1),par.time_scale,bsxfun(@minus,output.pphi_an(ind,:),output.pphi_an(ind,1))','k','displayname','$\sum_{j=R,Z,\varphi} \left(v_j^n \frac{\partial A_j^n}{\partial\varphi} \right)$');
            end
            
            % Left mid
            h(end+1)=plot(ha(2),par.time_scale,output.mm(ind,:),'b');
            
            % Left bottom
            h(end+1)=plot(ha(3),par.time_scale,output.Ekin(ind,:)/output.Ekin(ind,1)-1,'b');
            
            % Right (orbit)
            h(end+1)=plot(ha(4),squeeze(output.x   (ind,1,:)),squeeze(output.x   (ind,2,:)),'b','displayname','full orbit');
            h(end+1)=plot(ha(4),input.x(ind,1),input.x(ind,2),'rx','displayname','start point');
            %             h(5)=plot(ha(4),squeeze(output.x_gc(ind,1,:)),squeeze(output.x_gc(ind,2,:)),'r','displayname','gc orbit');
            
            % Title
            if exist('prec','var')
                title_1 = {['index = ',num2str(ind)],[' $\left<q\right>$ = ',num2str(prec.q_avg(ind))]};
                title_2= ['$\omega_b/\omega_d=$ ',num2str(prec.wb(ind)./prec.wd(ind),'%2.1e')];
            else
                title_1 = ['index = ',num2str(ind)];
                title_2 = '';
            end
            title(ha(4),title_1,'FontSize',20)
            title(ha(1),title_2,'FontSize',20)
            
            % Draw / delete
            drawnow;
            pause(0.2)
            
            h(h==0)=[]; %#ok<SAGROW>
            delete(h); clear('h')
        end
    case 2
        %% Phase space losses
        %% Distinquish particles
        expr_ind=~ejected & prec.pop.CO_PASSING;
%         expr_ind=true(size(ejected)) ;
        ind=find(expr_ind);
        
        plot_value=output.Delta_pphi(ind);
        %% Bins definition
        horizontal=prec.q_avg;
%         vertical=sign(vpll).*Eparralel;
        vertical=input.mm*B0./input.Ekin;
        
        x=horizontal(ind);
        y=vertical(ind);
        
        nr_x=100;
        x_0=min(x);
        x_1=max(x);
        delta_x=(x_1-x_0)/nr_x;
        
        xi=linspace(x_0,x_1,nr_x);
        Xedges = linspace(min(x)-0.5*delta_x,max(x)+0.5*delta_x,nr_x+1);
        
        nr_y=100;
        y_0=min(y); % Find minimum energy
        y_1=max(y); % Find maximum energy
        delta_y=(y_1-y_0)/nr_y;
        
        yi=linspace(y_0,y_1,nr_y);
        Yedges =linspace(y_0-0.5*delta_y,y_1+0.5*delta_y,nr_y+1); % Edges
        
        %% Binning
        % Binning        
        xr=interp1(Xedges,1:numel(Xedges),x,'linear'); % Index of edge (or try xi and nearest)
        yr=interp1(Yedges,1:numel(Yedges),y,'linear'); % Index of edge (or try yi and nearest)
        xr=floor(xr);
        yr=floor(yr);
        
        expr_fin=isfinite(xr) & isfinite(yr) & isfinite(plot_value);
        xr(~expr_fin)=1;
        yr(~expr_fin)=1;
        
        Z           =accumarray([yr xr],expr_fin,[nr_y nr_x]);                          % number of particles
        C           =accumarray([yr xr],expr_fin.*plot_value,[nr_y nr_x])./Z;      % Plot color averaged in bin
        
        % Correct scaling of color
        discr_pl_edges=linspace(-0.05,0.05,50);
        discr_pl_value=interp1(discr_pl_edges,1:numel(discr_pl_edges),plot_value,'nearest');
                
        C_total=NaN(size(C));
        for i=1:length(discr_pl_edges)
            expr_range= discr_pl_value==i;
            expr=expr_fin & expr_range;
            Z_i     =accumarray([yr xr],expr,[nr_y nr_x]);            
            C_total(Z_i>2)=sign(discr_pl_edges(i))*max(abs(plot_value(expr)));
        end
        
        EJ          =accumarray([yr xr],expr_fin & prec.ejected(ind),[nr_y nr_x]);      % ejected in 2D
        EJ_3D       =accumarray([yr xr],expr_fin & ejected(ind),[nr_y nr_x]);           % ejected in 3D
        TRAPPED     =accumarray([yr xr],expr_fin & prec.pop.ALL_TRAPPED(ind),[nr_y nr_x]); % trapped
        
        % Post process on matrices
        if any(Z(:)==0)
            warning('Not every bin has a particle')
        end
        
        % Ejected and Trapped as percentage
        EJ=EJ./Z*100;
        EJ_3D=EJ_3D./Z*100;
        TRAPPED=TRAPPED./Z*100;
        
        EJ_ALL=EJ==Z;
        EJ_ALL(Z==0)=false;
        
        %% Plotting
        hf=figure('tag',mfilename);
        ha=axes('parent',hf); hold(ha,'on')
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        xlabel(ha,'$\left<q\right>$','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        title_used=title_st;
        %         title_used{1}=[title_used{1},'\qquad %\varphi=$ ',num2str(mean(input.x(:,3))/pi,'%2.1f'),' $\pi$'];
        %         title(ha,{['sim \# ',par.ID,'\qquad $\varphi=$ ',num2str(mean(input.x(:,3))/pi,'%2.1f'),' $\pi$'],['$\Delta t=',num2str(par.dt,'%1.1e'),'$ [s], $t_{sim}=',num2str(par.time_scale(end),'%1.0e'),'$ [s]']},'FontSize',26)
        title(ha,title_used,'FontSize',26)
        set(ha,'FontSize',28);
        
%         [h]=imagescnan(xi,yi,C_total,'NaNColor',[1 1 1],'NaNMask',Z==0,'parent',ha);
        [~,h]=contourf(xi,yi,C,50,'linestyle','none','parent',ha);
%         [h]=imagesc(xi,yi,TRAPPED,'parent',ha);
        %[h]=imagesc(xi,yi,C,'parent',ha);
        
        cl=colorbar('peer',ha);
%         title(cl,'$\left<\frac{dp_\varphi}{dt}\right>$','FontSize',26,'interpreter','latex')
        title(cl,'$\Delta p_\varphi$','FontSize',26,'interpreter','latex')
        %        title(cl,'$\left<\frac{d\mu}{dt}\right>$','FontSize',26,'interpreter','latex')
        %         title(cl,'Lossed [\%]','FontSize',26,'interpreter','latex')
        axis xy
        xlim([x_0 x_1]); ylim([y_0 y_1]);
        %         caxis([-10 10])
        
        drawnow;
        %             delete(h); delete(cl);
    case 3
        %% Scatter plot the resonance behavior of trapped particles
        hf=figure('tag',mfilename);
        ha=axes('parent',hf);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title_used=title_st;
        title_used{1}=[title_used{1},' trapped particles'];
        
        gamma=abs(prec.wb./prec.wd);
        expr=gamma<10 & output.Delta_pphi>0;
        scatter(ha,gamma(expr),output.Delta_pphi(expr),'b.')
        
        set(ha,'FontSize',28)
        xlabel(ha,'$|\omega_b/\omega_d|$','FontSize',26)
        ylabel(ha,'$\Delta p_\varphi$','FontSize',26)
        title(ha,title_used,'FontSize',26)
        
    case 4
        %% Histogram of the resonance behavior of trapped particles (normalized)
        % Set gamma on x and Delta pphi on y
        gamma=(output.wb./output.wd);
        % Expression which are plotted
        expr=abs(gamma)<20;
        x=abs(gamma(expr));
        y=output.Delta_pphi(expr);
        
        nr_bins=[500 500];
        
        if verLessThan('matlab','8.6')
            % Plot using histndim with D and X - matrix:
            X=[x , y ];
            D(1,:)=nr_bins;
%             D(2,:)=[min(x) min(y)];
%             D(3,:)=[max(x) max(y)];
            D(2,:)=[0 -0.1];
            D(3,:)=[20 0.1];

            [Hor,Ver,Z ,N]=histndim(X,D);
            Z=Z(3:2:end-1,3:2:end-1);
            if any(Z'~=N)
                error('Some mistake in histndim?')
            end
            
            Xedges=Hor(2:2:end);
            Yedges=Ver(2:2:end);
            XC=0.5*(Xedges(2:end)+Xedges(1:end-1));
            YC=0.5*(Yedges(2:end)+Yedges(1:end-1));
            warning('Making use of histndim (ML File Exchange), which is not perfect in binning. Use a more recent version of Matlab for histcount2-function')
        else
            Xedges=linspace(0,20,nr_bins(1));
            Yedges=linspace(-0.1,0.1,nr_bins(2));
            [N,Xedges,Yedges] =     histcounts2(x,y,Xedges,Yedges); % Binning to check if this is correct
            XC=0.5*(Xedges(2:end)+Xedges(1:end-1));
            YC=0.5*(Yedges(2:end)+Yedges(1:end-1));
        end
        
        % Making a figure
        hf=figure('tag',mfilename);
        ha=axes('parent',hf);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title_used=title_st;
        title_used{1}=[title_used{1},' trapped particles'];
        
        % normalize the counts on the gamma-row number.
        Z=bsxfun(@times,N,100./sum(N,2))'; % Normalized on particle in same vertical box
        Z(Z<1e-1)=NaN;
        C=log10(Z); C(C==-Inf)=min(C(isfinite(C(:))));
        [~,h]=contourf(XC,YC,C,50,'linestyle','none','parent',ha);
                
        % Put in colorbar and labels
        xlabel(ha,'$|\omega_b/\omega_d|$','FontSize',26)
        ylabel(ha,'$\Delta p_\varphi/(Z_i e \psi_0)$','FontSize',26)
        zlabel(ha,'distribution in $|\omega_b/\omega_d|$ [\%]','FontSize',26)
        set(ha,'FontSize',28);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        cl=colorbar('peer',ha);
        title(cl,'\%','interpreter','latex','FontSize',26);
%         title(ha,title_used,'FontSize',26)
        set(h,'LineStyle','none');
        view(ha,2)
        y_lim=get(ha,'YLim');
        max_y_lim=0.1;%max(-y_lim(1),y_lim(2));
        set(ha,'YLim',[-1 1]*max_y_lim);
        set(ha,'XLim',[0 20]);
        % Alter ticks to have proper value (10^)
        order_high=2;
        order_low=floor(min(C(:)));
        if ~verLessThan('matlab','8.5')
            new_ticks=order_low:order_high;             % Values    of tick (e.g. -2 -1 0 1 2)
            new_tick_labels{length(new_ticks)}=[];      % Writing   of tick (e.g. 0.01 0.1 0 10 100) (Pre-allocate)
            for i=1:length(new_ticks)
                new_tick_labels{i}=['10^{',num2str(new_ticks(i)),'}'];
            end
            set(cl,'Ticks',new_ticks,'TickLabels',new_tick_labels)
        else
            title(cl,'\% 10\^{}','interpreter','latex','FontSize',26);
        end
        set(ha,'CLim',[-1 order_high])
        set(ha,'FontSize',28);
        
        expr_sum= XC>1.5 & XC<2.5;
        disp(['Sum as function of total N in |wb/wd|=2-peak ',num2str(sum(sum(N(expr_sum,:),2),1)./sum(N(:)))])
    case 5
        %% Figures in mu - pphi space with pphi the color (expr=dpphi_dt>0.01, color up to 0.04)
        expr= ~prec.ejected;
        %         delta_mm=abs(output.mm(:,end)-output.mm(:,1))./input.Ekin;
        %         expr=delta_mm>0.005 & ~prec.ejected;

        % Set the scaling of color
        %
        sc_data=output.Delta_pphi;
        sc_factor=(sc_data-0.00)/(0.04-0.00);
        %
%         sc_data=prec.wb./prec.wd;
%         sc_data(sc_data>0)=0; sc_data(sc_data<-15)=-15;
%         sc_factor=(sc_data)/(-10);
        %
        
        % Correct scaling of color
        sc_factor(sc_factor>1)=1;
        sc_factor(sc_factor<0)=0;
        sc_factor(~isfinite(sc_factor))=0;
        off_set_color=0;
        sc_factor=off_set_color+(1-off_set_color)*0.1*round(sc_factor*10);           % Use ten colors
        colors=unique(sc_factor(expr));
        
        horizontal_ini=input.pphi_kin(:,1);
        vertical_ini=input.mm(:,1)./input.Ekin*B0;
        horizontal_final=output.pphi_kin(:,end);
        vertical_final=output.mm(:,end)./input.Ekin*B0;
        
        % Figure trapped
        expr_type= expr & prec.pop.ALL_TRAPPED;
        hf=figure('tag',mfilename);
        ha=subplot(1,1,1,'parent',hf); hold(ha,'on')
        title_used=title_st;
        %
        xlabel(ha,'$p_\varphi/\psi_0$ ','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        set(ha,'FontSize',26);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title(ha,title_used,'FontSize',26)
        
        % Trapped
        expr_n=expr_type & ~ejected;
        for i=1:length(colors)
            expr_nn=expr_n & sc_factor==colors(i);
            if ~any(expr_nn); continue; end;
            h=plot(ha,horizontal_ini(expr_nn),vertical_ini(expr_nn),'g.','displayname','trapped');
            hasbehavior(h,'legend',false)
            C=get(h,'Color');
            set(h,'Color',C*colors(i));
        end
        plot(ha,NaN,NaN,'Color','g','Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        % Ejected
        expr_n=expr_type & ejected;
        for i=length(colors):-1:1
            expr_nn=expr_n & sc_factor==colors(i);
            if ~any(expr_nn); continue; end;
            h=plot(ha,horizontal_ini(expr_nn),vertical_ini(expr_nn),'y.','displayname','ejected');
            hasbehavior(h,'legend',false)
            C=get(h,'Color');
            set(h,'Color',C*colors(i));
        end
        plot(ha,NaN,NaN,'Color','y','Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        % Vectors
        expr_arrows=expr_type & ~ejected & output.Delta_pphi>0.01;
        nr_arrows=sum(expr_arrows);
        
        % Determine p1, where in the phase space they've gone
        p0=[horizontal_ini(expr_arrows)  ,   vertical_ini(expr_arrows)];
        p1=[horizontal_final(expr_arrows),   vertical_final(expr_arrows)];
        for i=1:100:nr_arrows
            h=vectarrow(p0(i,:),p1(i,:),ha,[0.494 0.184 0.556],'c');
            for j=1:length(h)
                hasbehavior(h(j),'legend',false);
            end
        end
        
        %% Figure CO-PASSING
        expr_type= expr & prec.pop.CO_PASSING;
        hf=figure('tag',mfilename);
        ha=axes('parent',hf); hold(ha,'on')
        title_used=title_st;
        
        xlabel(ha,'$p_\varphi/\psi_0$ ','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        set(ha,'FontSize',26);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title(ha,title_used,'FontSize',26)
        
        % Co-passing
        expr_n=expr & prec.pop.CO_PASSING & ~ejected;
        for i=1:length(colors)
            expr_nn=expr_n & sc_factor==colors(i);
            if ~any(expr_nn); continue; end;
            h=plot(ha,horizontal_ini(expr_nn),vertical_ini(expr_nn),'r.','displayname','co-passing');
            if i~=length(colors)
                hasbehavior(h,'legend',false)
            end
            C=get(h,'Color');
            set(h,'Color',C*colors(i));
        end
        
        % Ejected
        expr_n=expr_type & ejected;
        plot(ha,horizontal_ini(expr_n),vertical_ini(expr_n),'y.','displayname','ejected');
        
        % Vectors
        expr_arrows=expr_type & ~ejected & output.Delta_pphi>0.01;
        nr_arrows=sum(expr_arrows);
        
        % Determine p1, where in the phase space they've gone
        p0=[horizontal_ini(expr_arrows)  ,   vertical_ini(expr_arrows)];
        p1=[horizontal_final(expr_arrows),   vertical_final(expr_arrows)];
        for i=1:50:nr_arrows
            h=vectarrow(p0(i,:),p1(i,:),ha,[0.494 0.184 0.556],'c');
            for j=1:length(h)
                hasbehavior(h(j),'legend',false);
            end
        end
        
        %% Figure COUNTER-PASSING
        expr_type= expr & prec.pop.COUNTER_PASSING;
        hf=figure('tag',mfilename);
        ha=subplot(1,1,1,'parent',hf); hold(ha,'on')
        title_used=title_st;
        
        xlabel(ha,'$p_\varphi/\psi_0$ ','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        set(ha,'FontSize',28);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title(ha,title_used,'FontSize',26)
        
        % Counter-passing
        expr_n=expr & prec.pop.COUNTER_PASSING;
        for i=1:length(colors)
            expr_nn=expr_n & sc_factor==colors(i) & ~ejected;
            if ~any(expr_nn); continue; end;
            h=plot(ha,horizontal_ini(expr_nn),vertical_ini(expr_nn),'b.','displayname','counter-passing');
            if i~=length(colors)
                hasbehavior(h,'legend',false)
            end
            C=get(h,'Color');
            set(h,'Color',C*colors(i));
        end
        
        % Ejected
        expr_n=expr_type & ejected;
        h=plot(ha,horizontal_ini(expr_n),vertical_ini(expr_n),'y.','displayname','ejected');
        
        % Vectors
        expr_arrows=expr_type & ~ejected & output.Delta_pphi>0.01;
        nr_arrows=sum(expr_arrows);
        
        % Determine p1, where in the phase space they've gone
        p0=[horizontal_ini(expr_arrows)  ,   vertical_ini(expr_arrows)];
        p1=[horizontal_final(expr_arrows),   vertical_final(expr_arrows)];
        for i=1:50:nr_arrows
            h=vectarrow(p0(i,:),p1(i,:),ha,[0.494 0.184 0.556],'c');
            for j=1:length(h)
                hasbehavior(h(j),'legend',false);
            end
        end
        
    case 6
        %% Scatter plot in mu/E and <q> phase space
        hf=figure('tag',mfilename);
        ha=subplot(1,2,1,'parent',hf); hold(ha,'on')
        title_used=title_st;
        
        xlabel(ha,'$\left<q\right>$ (2D)','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        set(ha,'FontSize',26);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title(ha,title_used,'FontSize',26)
        
        % General parameters
        expr=true(size(ejected)) & (input.Ekin>40*1e3)&(input.Ekin<50*1e3);       % Which to plot
        expr_passing=expr & prec.pop.COUNTER_PASSING;
        horizontal=prec.q_avg;
        vertical=input.mm*B0./input.Ekin;
        
        % Scatter plots
        %          expr_n=expr & prec.ejected & (input.v(:,3)>0);
        %         h=plot(ha,horizontal(expr_n),vertical(expr_n),'k.','displayname','ejected in 2D');
        %         hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        % Plot trapped
        expr_n=expr & ~ejected & prec.pop.ALL_TRAPPED;
        h=plot(ha,horizontal(expr_n),vertical(expr_n),'g.','displayname','trapped and confined');
        set(h,'Color',[0 0.9 0]);
        hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        
        expr_n=expr & ejected & ~prec.ejected & prec.pop.ALL_TRAPPED;
        h=plot(ha,horizontal(expr_n),vertical(expr_n),'y.','displayname','trapped and ejected in 3D');
        set(h,'Color',[0.9 0.9 0]);
        hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        
        % Plot passing
        expr_n=expr_passing & ~ejected;
        h=plot(ha,horizontal(expr_n),vertical(expr_n),'b.','displayname','co-current and confined');
        set(h,'Color',[0 0 0.9]);
        hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        
        expr_n=expr_passing & ejected;
        h=plot(ha,horizontal(expr_n),vertical(expr_n),'m.','displayname','co-current and ejected in 3D');
        set(h,'Color',[0.494 0.184 0.556]);
        hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        
        
        % Vectors
        %         expr_n=expr_passing & ejected;
        %         nr_arrows=sum(expr_n);
        
        % Fit <q> to pphi to determine Delta <q>
        %         expr_fit=isfinite(output.pphi_kin(:,1)) & isfinite(prec.q_avg) & expr_passing;
        %         p=polyfit(output.pphi_kin(expr_fit,1),prec.q_avg(expr_fit),2);
        %         q0=polyval(p,output.pphi_kin(:,1));
        %         q1=polyval(p,output.pphi_kin(:,end));
        %         delta_q=q1-q0;
        
        %         p0=[horizontal(expr_n)  ,   vertical(expr_n)];
        %         p1=[horizontal(expr_n)+delta_q(expr_n)  ,   output.mm(expr_n,end)./input.Ekin(expr_n)]*1.02;
        
        %         p0=[output.pphi_kin(expr_n,1)  ,   vertical(expr_n)];
        %         p1=[output.pphi_kin(expr_n,end),   output.mm(expr_n,end)./input.Ekin(expr_n)]*1.02;
        %         for i=1:15:nr_arrows
        %             h=vectarrow(p0(i,:),p1(i,:),ha,[0.5 0.5 0.5],'r');
        %             for j=1:length(h)
        %                 hasbehavior(h(j),'legend',false);
        %             end
        %         end
        
        % Legend
        hl=legend(ha,'show');
        set(hl,'interpreter','latex','location','southwest')
    case 7
        %% Scatter plot in mu/E and <q> phase space
        hf=figure('tag',mfilename);
        ha=subplot(1,1,1,'parent',hf); hold(ha,'on')
        title_used=title_st;
        
        xlabel(ha,'$p_\varphi/\psi_0$','FontSize',26)
        ylabel(ha,'$\mu/E_{kin}$ [1/T]','FontSize',26)
        set(ha,'FontSize',26);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title(ha,title_used,'FontSize',26)
        
        % General parameters
        expr=~prec.ejected & any(isfinite(output.pphi_kin),2) & (input.Ekin>40*1e3)&(input.Ekin<50*1e3) & any(input.mm(:,:)~=0,2);       % Which to plot
        expr_passing = expr & prec.pop.COUNTER_PASSING;
        horizontal=output.pphi_kin(:,1);
        vertical=input.mm./input.Ekin;
        
        % Scatter plots
        expr_n=expr & prec.ejected & (input.v(:,3)<0);
        if any(expr_n)
            h=plot(ha,horizontal(expr_n),vertical(expr_n),'k.','displayname','ejected in 2D (with $v_\parallel<0$');
            hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        end
        expr_n=expr & ~ejected & prec.pop.ALL_TRAPPED;
        if any(expr_n)
            h=plot(ha,horizontal(expr_n),vertical(expr_n),'g.','displayname','trapped and confined');
            set(h,'color',[0 0.9 0]);
            hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        end
        expr_n=expr & ejected & ~prec.ejected & prec.pop.ALL_TRAPPED;
        if any(expr_n)
            h=plot(ha,horizontal(expr_n),vertical(expr_n),'y.','displayname','trapped and ejected in 3D');
            set(h,'color',[0.9 0.9 0]);
            hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        end
        expr_n=expr_passing & ~ejected ;
        if any(expr_n)
            h=plot(ha,horizontal(expr_n),vertical(expr_n),'b.','displayname','co-current and confined');
                    set(h,'color',[0 0 0.9]);
%             set(h,'color',[0.9 0 0]);
            hasbehavior(h,'legend',false); plot(NaN,NaN,'Color',get(h,'Color'),'Marker',get(h,'Marker'),'LineStyle','none','displayname',get(h,'displayname'),'MarkerSize',50); % Make fake point (marker) for legend
        end
        % Vectors
        expr_n=expr_passing & ejected;
        nr_arrows=sum(expr_n);

        part_indexes=find(expr_n);
        p0=[output.pphi_kin(part_indexes,1)  ,   vertical(part_indexes,1)];
        
        last=isfinite(output.pphi_kin(expr_n,:));
        for i=1:50:nr_arrows
            part_ind=part_indexes(i);
            last_ind=find(last(i,:),1,'last');
            
            p1=[output.pphi_kin(part_ind,last_ind),    output.mm(part_ind,last_ind)./input.Ekin(part_ind)];
            h=vectarrow(p0(i,:),p1,ha,[0.494 0.184 0.556],'c');
            for j=1:length(h)
                hasbehavior(h(j),'legend',false);
            end
        end
        
        % Legend
        hl=legend(ha,'show');
        set(hl,'interpreter','latex','location','southwest')
    case 8
        %% Scatter plot in mu/E and <q> phase space
        hf=figure('tag',mfilename);
        ha=subplot(1,1,1,'parent',hf); hold(ha,'on')
        xlabel(ha,'$p_\varphi/(Z_i e \psi_0)$','FontSize',26)
        ylabel(ha,'$\mu B_0/E_{kin}$','FontSize',26)
        set(ha,'FontSize',28);
        if ~verLessThan('matlab','8.5');	set(ha,'TickLabelInterpreter','latex');	end
        title_used=title_st;
        min_particles_box=10;
        %title_used{1}=[title_used{1},', lossses with minimal ',num2str(min_particles_box),' per bin'];
        %title(ha,title_used,'FontSize',26)
        
        % General parameters
%         expr=true(input.N_total,1);       % Which to plot
        expr=prec.pop.ALL_TRAPPED | prec.pop.ALL_PASSING;
        expr_passing = expr & prec.pop.ALL_PASSING;
        x=input.pphi_kin(:,1);
        y=input.mm(:,1)*B0./input.Ekin;
        
        nr_bins=[250 250];
        
        x(~expr,:)=NaN;
        y(~expr,:)=NaN;
        if verLessThan('matlab','8.6')
            X=[x , y ];
            D(1,:)=nr_bins;
%             xl=[min(x),max(x)];
%             yl=[min(y) max(y)];
            xl=[-1.320557,0.21213];
            yl=[0 1.26];

            D(2,:)=[xl(1) yl(1)];
            D(3,:)=[xl(2) yl(2)];
            [Hor,Ver,Z ,N]=histndim(X,D);
            Z=Z(3:2:end-1,3:2:end-1);
            if any(Z'~=N)
                error('Some mistake in histndim?')
            end
            
            XC=Hor(2:2:end-1);
            YC=Ver(2:2:end-1);
            
            % PDF 
            ZPDF=Z/(mean(diff(XC))*mean(diff(YC))*sum(expr));
            
            % Ejected
            [~,~,~,NE]=histndim(X(ejected,:),D);
            ZE=NE';
            
             % Ripple trapped particles
            expr_ripple_trapped=output.nr_vpll_crossing>3 & output.nr_midplane_crossing==0;
            [~,~,~,Nripple]=histndim(X(expr_ripple_trapped,:),D);
            Zripple=Nripple';
            
            % Synergy
            if exist('ej_AB','var') && exist('ej_ST','var') && exist('ej_RMP','var')
                expr_s_plus= ej_AB & ~ej_ST & ~ej_RMP;
                expr_s_min =~ej_AB &  ej_ST &  ej_RMP;
                
                [~,~,~,N_s_plus]=histndim(X(expr_s_plus,:),D);
                Z_s_plus=N_s_plus';       
                [~,~,~,N_s_min]=histndim(X(expr_s_min,:),D);
                Z_s_min=N_s_min';        
                Z_s=Z_s_plus-Z_s_min;
            end
            
            % Trapped and Passing
            [~,~,~,NCOP]=histndim(X(prec.pop.CO_PASSING,:),D);          % ZCOP, number of co-passing particles
            [~,~,~,NCOUNTP]=histndim(X(prec.pop.COUNTER_PASSING,:),D);  % ZCOUNP, number of counter-passing particles
            [~,~,~,NTR]=histndim(X(prec.pop.ALL_TRAPPED,:),D);          % ZTR, number of trapped particles
            [~,~,~,NPA]=histndim(X(prec.pop.ALL_PASSING,:),D);          % ZPA, number of passing particles
            ZCOP=logical(NCOP)';
            ZCOUNP=logical(NCOUNTP)';
            ZTR=logical(NTR)';
            ZPA=logical(NPA)';
            
            warning('Making use of histndim (ML File Exchange), which is not perfect in binning. Use a more recent version of Matlab for histcount2-function')
        else    % Making use of histcount2
            % Normal count and  edges
            if exist('xl','var')
                Xedges=linspace(xl(1),xl(2),nr_bins(1)+1);
                Yedges=linspace(yl(1),yl(2),nr_bins(2)+1);
                [N,Xedges,Yedges] =     histcounts2(x,y,Xedges,Yedges); % Main binning algorithmparticles 
            else
                [N,Xedges,Yedges] =     histcounts2(x,y,nr_bins); % Main binning algorithmparticles 
            end
            Z=N';
            
            XC=linspace(Xedges(1),Xedges(end),nr_bins(1)*2+1);  % Vector of 
            XC=XC(2:2:end-1);
            YC=linspace(Yedges(1),Yedges(end),nr_bins(2)*2+1);
            YC=YC(2:2:end-1);
            % Probability density function
            [NPDF,Xedges,Yedges] =     histcounts2(x,y,Xedges,Yedges,'Normalization','pdf'); % 
            ZPDF=NPDF';
            
            % Ejected particles
            [NE] =     histcounts2(x(ejected,:),y(ejected,:),Xedges,Yedges); 
            ZE=NE';
            
             % Ripple trapped particles
%             expr_ripple_trapped=output.nr_vpll_crossing>3 & output.nr_midplane_crossing==0;
%             [Nripple] =     histcounts2(x(expr_ripple_trapped,:),y(expr_ripple_trapped,:),Xedges,Yedges); 
%             Zripple=Nripple';
            
            % Synergy
            if exist('ej_AB','var') && exist('ej_ST','var') && exist('ej_RMP','var')
                expr_s_plus= ej_AB & ~ej_ST & ~ej_RMP;
                expr_s_min =~ej_AB &  ej_ST &  ej_RMP;
                [N_s_plus] =     histcounts2(x(expr_s_plus,:),y(expr_s_plus,:),Xedges,Yedges); 
                Z_s_plus=N_s_plus';        
                [N_s_min] =     histcounts2(x(expr_s_min,:),y(expr_s_min,:),Xedges,Yedges); 
                Z_s_min=N_s_min';        
                Z_s=Z_s_plus-Z_s_min;
            end
            
            % Trapped and Passing
            [NCOP] =    histcounts2(x(prec.pop.CO_PASSING,:)        ,y(prec.pop.CO_PASSING,:)       ,Xedges,Yedges); % ZCOP, number of co-passing particles
            [NCOUNTP] = histcounts2(x(prec.pop.COUNTER_PASSING,:)   ,y(prec.pop.COUNTER_PASSING,:)  ,Xedges,Yedges); % ZCOUNP, number of counter-passing particles
            [NTR] =     histcounts2(x(prec.pop.ALL_TRAPPED,:)       ,y(prec.pop.ALL_TRAPPED,:)      ,Xedges,Yedges); % ZTR, number of trapped particles
            [NPA] =     histcounts2(x(prec.pop.ALL_PASSING,:)       ,y(prec.pop.ALL_PASSING,:)      ,Xedges,Yedges); % ZPA, number of passing particles
            ZCOP=logical(NCOP)';
            ZCOUNP=logical(NCOUNTP)';
            ZTR=logical(NTR)';  % Convert to logical for boundary
            ZPA=logical(NPA)';  % Convert to logical for boundary
        end
        
        % Contraint of particles in box!
        NaN_expr=Z<min_particles_box | Z==0;
        Z(NaN_expr)=NaN; 
        ZE(NaN_expr)=NaN; 
%         Zripple(NaN_expr)=NaN; 
        ZPDF(NaN_expr)=NaN; 
        if exist('Z_s_plus','var') || exist('Z_s_min','var') 
            Z_s_plus(NaN_expr)=NaN;
            Z_s_min (NaN_expr)=NaN;
            Z_s     (NaN_expr)=NaN;
            Z_s_plus=Z_s_plus./Z*100;
            Z_s_min =Z_s_min ./Z*100;
            Z_s     =Z_s     ./Z*100;
        end
        % Calculate fraction of ejected particles
        ZE=ZE./Z*100;
%         Zripple=Zripple./Z*100;
        
        % Make contour plot
%         C=log10(ZE); C(C==-Inf)=min(C(isfinite(C(:))));
        C=Z_s;
        [~,hc]=contourf(XC,YC,C,100,'parent',ha,'linestyle','none');
%         [hc,hcnan]=imagescnan(XC,YC,C,'parent',ha); set(ha,'YDir','normal');
%         set(hc,'displayname','losses as percentage')
        hasbehavior(hc   ,'legend',false);
        if exist('hcnan','var')
            hasbehavior(hcnan,'legend',false);
        end
        caxis(ha,'auto'); 
        
        %  TRAPPED CONTOUR
        [hg]=contourc(XC,YC,double(ZTR),1);
        i=1;
        while i~=size(hg,2)+1
            ind_advance=hg(2,i);
            hl=plot(ha,hg(1,i+1:i+ind_advance),hg(2,i+1:i+ind_advance),'color','g','linewidth',1,'displayname','Area of trapped particles');
            if i~=1
                hasbehavior(hl,'legend',false);
            end
            i=i+ind_advance+1;
            % Close contour if only 1 contour
            if ind_advance==size(hg,2)-1
                X=get(hl,'XData'); Y=get(hl,'YData');
                set(hl,'XData',[X,X(1)],'YData',[Y,Y(1)])
            end
        end
%         % PASSING CONTOUR
%         [hr]=contourc(XC,YC,double(ZPA),1);
%         i=1;
%         while i~=size(hr,2)+1
%             ind_advance=hr(2,i);
%             hl=plot(ha,hr(1,i+1:i+ind_advance),hr(2,i+1:i+ind_advance),'color','r','linewidth',1,'displayname','Area of passing particles');
%             if i~=1
%                 hasbehavior(hl,'legend',false);
%             end
%             i=i+ind_advance+1;
%         end        
        
       %  CO-PASSING CONTOUR
        [hr]=contourc(XC,YC,double(ZCOP),1);
        i=1;
        while i~=size(hr,2)+1
            ind_advance=hr(2,i);
            hl=plot(ha,hr(1,i+1:i+ind_advance),hr(2,i+1:i+ind_advance),'color','k','linewidth',1,'displayname','Area of counter-current passing particles');
            if i~=1
                hasbehavior(hl,'legend',false);
            end
            i=i+ind_advance+1;
        end        
%         
        % COUNTER PASSING CONTOUR
        [hr]=contourc(XC,YC,double(ZCOUNP),1);
        i=1;
        while i~=size(hr,2)+1
            ind_advance=hr(2,i);
            hl=plot(ha,hr(1,i+1:i+ind_advance),hr(2,i+1:i+ind_advance),'color','r','linewidth',1,'displayname','Area of co-current passing particles');
            if i~=1
                hasbehavior(hl,'legend',false);
            end
            i=i+ind_advance+1;
        end        
        
        
        % Put in colorbar and labels
        cl=colorbar('peer',ha);
        title(cl,'\%','interpreter','latex','FontSize',26);
        
        % Make colordata on logarithmic scale
        % Alter ticks to have proper value (10^)
%         order_high=2;
%         order_high=ceil(max(C(:)));
%         order_low=floor(min(C(:)));
%         if ~verLessThan('matlab','8.5')
%             new_ticks=order_low:order_high;             % Values    of tick (e.g. -2 -1 0 1 2)
%             new_tick_labels{length(new_ticks)}=[];      % Writing   of tick (e.g. 0.01 0.1 0 10 100) (Pre-allocate)
%             for i=1:length(new_ticks)
%                 % new_tick_labels{i}=num2str(10^new_ticks(i),'%3.3f');%['10^{',num2str(new_ticks(i)),'}'];
%                 new_tick_labels{i}=['10^{',num2str(new_ticks(i)),'}'];
%             end
%             set(cl,'Ticks',new_ticks,'TickLabels',new_tick_labels)
%         else
%             title(cl,'\% 10\^{}','interpreter','latex','FontSize',26);
%         end
% %         set(ha,'Clim',[0 order_high])
% %         set(cl,'FontSize',28);
        xlim(ha,[min(XC),max(XC)])
        ylim(ha,[min(YC),max(YC)])
        % Legend
   %     hl=legend(ha,'show');
%         set(hl,'interpreter','latex','location','southwest')
end