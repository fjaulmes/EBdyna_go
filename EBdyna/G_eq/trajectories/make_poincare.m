
%% MAKE A FIGURE FOR POINCARE-PLOT
function [ha,hq,psi_to_psi_overline]=make_poincare
global par dim maps
q_to_psi = @(q) interp1(1:length(dim.psi_scale),dim.psi_scale,interp1(dim.q_initial_profile,1:length(dim.q_initial_profile),q),'cubic'); % psi from q
psi_to_q = @(psi) interp1(dim.psi_scale,dim.q_initial_profile,psi,'cubic');
psi_to_psi_overline=[];
% Plot lines of rational values with RMP_mode
RMP_mode=2;
rational_q_surfaces=4:(-1/RMP_mode):1; % For q=1 to q=4

%% Make a figure
try set(groot,'defaulttextinterpreter','latex');    catch err; end;
delete(findall(0,'type','figure','tag','Poincare_map'))
hf=figure('Name','Poincare map','tag','Poincare_map');

open_fig_ST_xy='./2016-12-13_ST_xy_frame_4.fig';
open_fig_ST_RZ='./2016-12-13_ST_RZ_frame_4.fig';

switch par.poincare_type
    case 'sfl'
        %% psi-theta=space
        psi_to_psi_overline = @(psi) 1-psi/maps(1).psi_global;   % normalized psi from psi
        q_to_psi_overline = @(q) psi_to_psi_overline(q_to_psi(q));                  % normalized psi from q
        
        [AX,h1,h2]=plotyy([0 1],[0 1],[0 1],[0 1],'parent',hf);                     % Make an double y-axis figure
        delete(h1); delete(h2); linkaxes(AX,'xy'); clear h1 h2
        ha=AX(1);   hold(ha,'on'); set(ha,'FontSize',20);                           % Have an psi-theta plot figure
        hq=AX(2);   hold(hq,'on'); set(hq,'FontSize',20);                           % Have an q-theta plot figure
        set(ha,'YColor','k'); set(hq,'YColor','r');
        xlabel(ha,'$\theta$ [rad]','FontSize',24,'interpreter','latex');
        ylabel(ha,'$\overline\psi$','FontSize',24,'interpreter','latex');
        ylabel(hq,'$q$','FontSize',24,'interpreter','latex');
        set(ha,'XLim',[0 2*pi],'Ylim',[0 psi_to_psi_overline(0)]);
        
        % Determine psi on these surfaces
        psi_overline_q_surfaces=q_to_psi_overline(rational_q_surfaces);
        
        % Make lines on q-values and label them
        q_surf_name{length(rational_q_surfaces)}=[];
        for i=1:length(rational_q_surfaces)
            if mod(rational_q_surfaces(i),1)==0
                line([0 2*pi],[psi_overline_q_surfaces(i) psi_overline_q_surfaces(i)],'color','r','linestyle','-','LineWidth',1.5,'parent',hq);
                q_surf_name{i}=[num2str(rational_q_surfaces(i)*RMP_mode),'/',num2str(RMP_mode)];
            else
                line([0 2*pi],[psi_overline_q_surfaces(i) psi_overline_q_surfaces(i)],'color','r','linestyle',':','LineWidth',1.5,'parent',hq);
            end
        end
        set(hq,'YTick',fliplr(psi_overline_q_surfaces),'YtickLabel',fliplr(q_surf_name))
        set(ha,'YTick',fliplr(psi_overline_q_surfaces));
        
    case 'tor'
        hq=[];
        q_XZ=psi_to_q(maps(1).psi_XZ);
        
        ha=axes('parent',hf);
        hold(ha,'on'); set(ha,'FontSize',20);
        axis(ha,'equal','xy')
        xlabel(ha,'$R$','FontSize',24,'interpreter','latex');
        ylabel(ha,'$Z$','FontSize',24,'interpreter','latex');
        if par.APPLY_SAWTOOTH && exist(open_fig_ST_RZ,'file')
            hf_temp=openfig(open_fig_ST_RZ);
            ha_temp=findall(hf_temp,'type','axes');
            ha_temp=ha_temp(2);
            h_temp=findall(ha_temp,'type','contour');
            set(flipud(h_temp),'parent',ha);
            close(hf_temp)
            colormap(ha,'jet')
            caxis(ha,[0 max(h_temp.ZData(:))]);
        end
        % Plot LCFS
        c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps(1).psi_norm_XZ',[513 513]);
        plot(ha,c(1,2:end),(c(2,2:end)),'r','displayname','LCFS','LineWidth',2)
        contour(dim.R0+dim.scale_X,dim.scale_Z,q_XZ',rational_q_surfaces,'ShowText','on','parent',ha,'color','k','linewidth',1);
        %         for i=1:length(rational_q_surfaces)
        %             c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps(1).psi_XZ',[psi_q_surfaces(i) psi_q_surfaces(i)]);
        %             if mod(rational_q_surfaces(i),1)==0
        %                 arg_line_color_style='r';
        %             else
        %                 arg_line_color_style='r:';
        %             end
        %             plot(ha,c(1,2:end),(c(2,2:end)),arg_line_color_style,'displayname',['resonant surface, q= ',[num2str(rational_q_surfaces(i)*RMP_mode),'/',num2str(RMP_mode)]],'LineWidth',1)
        %         end
    case 'xy'
        hq=[];
        q_XZ=psi_to_q(maps(1).psi_XZ);
        
        ha=axes('parent',hf);
        hold(ha,'on'); set(ha,'FontSize',20);
        axis(ha,'equal','xy')
        xlabel(ha,'$x$','FontSize',24,'interpreter','latex');
        ylabel(ha,'$y$','FontSize',24,'interpreter','latex');
        if par.APPLY_SAWTOOTH && exist(open_fig_ST_RZ,'file')
            hf_temp=openfig(open_fig_ST_xy);
            ha_temp=findall(hf_temp,'type','axes');
            ha_temp=ha_temp(2);
            h_temp=findall(ha_temp,'type','contour');
            %             set(flipud(h_temp),'parent',ha);
            [~,h_new]=contourf(h_temp.XData,h_temp.YData,h_temp.ZData,200,'linestyle',':','linecolor','r','linewidth',1,'parent',ha);
            close(hf_temp)
            colormap(ha,'summer')
            caxis(ha,[0 max(h_new.ZData(:))]);
            xlim([-1 1]*max(abs(h_new.XData(:))));
            ylim([-1 1]*max(abs(h_new.YData(:))));
        end
        % Plot q-values
        r=sqrt(2*maps(1).chi_XZ);
        theta= maps(1).theta_normal_XZ;
        x=r.*cos(theta);
        y=r.*sin(theta);
        contour(x,y,q_XZ,rational_q_surfaces,'ShowText','on','parent',ha,'color','b','linewidth',1);
    otherwise
        error('plotting of poincare not specified')
end

end
