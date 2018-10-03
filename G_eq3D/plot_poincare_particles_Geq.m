clear all
clc
format compact
format long
initialize_folder_names;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename,'scale_X','scale_Z','R0');
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename,'psi_scale','rho_tor_scale','psi_global');
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename,'q_initial_profile');
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename,'theta_XZsmall_map','psi_XZsmall_map','q_initial_XZsmall_map');

FILENUMBER=1
SUBNUMBER=1

psi_limit1=interp1(rho_tor_scale,psi_scale,0.7)
psi_limit2=interp1(rho_tor_scale,psi_scale,0.9)

ID='1915577'

figure;
hold on
contour(scale_X+R0,scale_Z,psi_XZsmall_map',[0 0],'y','linewidth',2);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[1.5 1.5],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[2.0 2.0],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[2.5 2.5],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[3.0 3.0],'r','linewidth',1);
axis xy equal square
theta_ej=[];
psi_ej=[];
for file_index=1:SUBNUMBER:FILENUMBER
    filename=['output/G_eq_19xx_full_process' num2str(file_index) '.mat']
    try
    load(filename,'output','ejected')
    part_count=0;
    
    for part_index=1:size(output.x_gc,1)
        valid_ts=find(~isnan(output.x_gc(part_index,1,:)));
        if ejected(part_index)
            colorstring='m.';
            theta_ej=[theta_ej ...
                squeeze(interp2(scale_X+R0,scale_Z,theta_XZsmall_map',output.x_gc(part_index,1,valid_ts),output.x_gc(part_index,2,valid_ts),'*linear'))'];
            psi_ej=[psi_ej ...
                squeeze(interp2(scale_X+R0,scale_Z,psi_XZsmall_map',output.x_gc(part_index,1,valid_ts),output.x_gc(part_index,2,valid_ts),'*linear'))'];
        else
            colorstring='k.';
        end
        part_count=part_count+length(valid_ts);
        psi_gc=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',output.x_gc(part_index,1,1),output.x_gc(part_index,2,1),'*linear');
        if part_index<round(0.25*size(output.x_gc,1))
            Rgc=squeeze(output.x_gc(part_index,1,valid_ts));
            Zgc=squeeze(output.x_gc(part_index,2,valid_ts));
            plot(Rgc,Zgc,colorstring);
        elseif part_index<round(0.5*size(output.x_gc,1))
            if psi_gc>psi_limit1
                psi_gc
                Rgc=squeeze(output.x_gc(part_index,1,valid_ts));
                Zgc=squeeze(output.x_gc(part_index,2,valid_ts));
                plot(Rgc,Zgc,colorstring);
            end
        else
            
            if psi_gc>psi_limit2
                psi_gc
%                 for ts=1:SUBTIME:length(valid_ts)
%                     plot(output.x_gc(part_index,1,valid_ts(ts)),output.x_gc(part_index,2,valid_ts(ts)),'k.');
%                 end
                Rgc=squeeze(output.x_gc(part_index,1,valid_ts));
                Zgc=squeeze(output.x_gc(part_index,2,valid_ts));
                plot(Rgc,Zgc,colorstring);
            end
            
        end
    end
    catch
        disp('COULD NOT LOAD FILE!')
    end
end

contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[1.5 1.5],'r--','linewidth',1.1);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[2.0 2.0],'r','linewidth',1.3);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[2.5 2.5],'r--','linewidth',1.1);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[3.0 3.0],'r','linewidth',1.3);
psi_r95=interp1(rho_tor_scale,psi_scale,0.95)
contour(scale_X+R0,scale_Z,psi_XZsmall_map',[psi_r95 psi_r95],'b','linewidth',1.5);

disp('move to radial, poloidal plot ....')

%%
psi_norm_scale=1+psi_scale/psi_global

figure;
hold on

ej_count=0;
for file_index=1:SUBNUMBER:FILENUMBER
%     if file_index~=73
%     filename=['output_' num2str(ID) '/G_eq_19xx_full_process' num2str(file_index) '.mat']
%     end
    load(filename,'output','ejected')
    ej_count=ej_count+length(find(ejected));
    part_count=0;
    output.theta=reshape(output.theta,size(output.theta,1)*size(output.theta,2),1);
    output.theta=output.theta(~isnan(output.theta));
    output.psi=reshape(output.psi,size(output.psi,1)*size(output.psi,2),1);
    output.psi=output.psi(~isnan(output.psi));
    output_psi_bar=interp1(psi_scale,psi_norm_scale,output.psi);
    plot(output.theta,output_psi_bar,'k.');
    
%     for part_index=1:size(output.x_gc,1)
%         ini_index=part_count+1;
%         valid_ts=find(~isnan(output.x_gc(part_index,1,:)));
%         part_count=part_count+length(valid_ts);
%         end_index=part_count;

%         psi_gc=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',output.x_gc(part_index,1,1),output.x_gc(part_index,2,1),'*linear');
%         if part_index<round(0.25*size(output.x_gc,1))
%             for pindex=ini_index:end_index
%                 plot(output.theta(pindex),output.psi(pindex),'k.');
%             end
%         elseif part_index<round(0.5*size(output.x_gc,1))
%             if psi_gc>psi_limit1
%                 psi_gc
%                 for pindex=ini_index:end_index
%                     plot(output.theta(pindex),output.psi(pindex),'k.');
%                 end
%             end
%         else
%             
%             if psi_gc>psi_limit2
%                 psi_gc
%                 for pindex=ini_index:end_index
%                     plot(output.theta(pindex),output.psi(pindex),'k.');
%                 end
%             end
%             
% %         end
%     end
end
ej_count
xlabel('\theta')
ylabel('$\bar{\psi}$','interpreter','latex')


%%
output_psi_ej_bar=interp1(psi_scale,psi_norm_scale,psi_ej);

plot(theta_ej,output_psi_ej_bar,'m.');

    
psi_q1=interp1(q_initial_profile,psi_norm_scale,1.5)
psi_q2=interp1(q_initial_profile,psi_norm_scale,2)
psi_q3=interp1(q_initial_profile,psi_norm_scale,2.5)
psi_q4=interp1(q_initial_profile,psi_norm_scale,3)
psi_q5=interp1(q_initial_profile,psi_norm_scale,3.5)

plot([0 2*pi],[psi_q1 psi_q1],'r--','linewidth',1.3)
plot([0 2*pi],[psi_q2 psi_q2],'r','linewidth',1.5)
plot([0 2*pi],[psi_q3 psi_q3],'r--','linewidth',1.3)
plot([0 2*pi],[psi_q4 psi_q4],'r','linewidth',1.5)
plot([0 2*pi],[psi_q5 psi_q5],'r--','linewidth',1.3)

xlim([0 2*pi])
% ylim([-0.0109 0])

%%

psi_r88=interp1(rho_tor_scale,psi_norm_scale,0.88)
psi_r90=interp1(rho_tor_scale,psi_norm_scale,0.9)
psi_r95=interp1(rho_tor_scale,psi_norm_scale,0.95)

plot([0 2*pi],[psi_r95 psi_r95],'b','linewidth',2)

ylim([psi_r88 1])
set(gca,'fontsize',26)