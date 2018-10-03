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

FILENUMBER=256
SUBNUMBER=1

psi_limit1=interp1(rho_tor_scale,psi_scale,0.7)
psi_limit2=interp1(rho_tor_scale,psi_scale,0.9)

ID='1915611'

%%
psi_norm_scale=1+psi_scale/psi_global;

figure;
hold on
DL_values=[];
ej_count=0;
DLMAX=0;
theta_ej=[];
psi_ej=[];
try
    load DL_values.mat
catch
    for file_index=1:SUBNUMBER:FILENUMBER
        try
            filename=['./output_' num2str(ID) '/G_eq_19xx_full_process' num2str(file_index) '.mat'];
            load(filename,'output','ejected')
            DL_values=[DL_values;output.Delta_l];
            DLMAX=max(DL_values);
        end
    end
    save DL_values.mat DLMAX DL_values
end
%%
cmhot=hot(256);
for file_index=1:SUBNUMBER:FILENUMBER
    try
        filename=['./output_' num2str(ID) '/G_eq_19xx_full_process' num2str(file_index) '.mat']
        load(filename,'output','ejected')
        ej_count=ej_count+length(find(ejected));
        for part_index=1:size(output.x_gc,1)
            valid_ts=find(~isnan(output.x_gc(part_index,1,:)));
            output.theta(part_index,valid_ts)=interp2(scale_X+R0,scale_Z,theta_XZsmall_map',output.x_gc(part_index,1,valid_ts),output.x_gc(part_index,2,valid_ts),'*linear');
            if ejected(part_index)
                %             colorstring='m.';
                colorindex=256-round(255*output.Delta_l(part_index)/DLMAX);
                colorstring=[cmhot(colorindex,1) cmhot(colorindex,2) cmhot(colorindex,3)];
                theta_ej=[theta_ej ...
                    squeeze(output.theta(part_index,valid_ts))];
                psi_ej=[psi_ej ...
                    squeeze(output.psi(part_index,valid_ts))];
            else
%                 colorindex=256-round(255*output.Delta_l(part_index)/DLMAX);
%                 colorstring=[cmhot(colorindex,1) cmhot(colorindex,2) cmhot(colorindex,3)]; 
                colorstring=[0 0 0];
            end
            output_psi_bar=interp1(psi_scale,psi_norm_scale,output.psi(part_index,valid_ts),'*PCHIP');
            plot(output.theta(part_index,valid_ts),output_psi_bar,'color',colorstring,'linestyle','none','marker','.');
        end
    end
end
ej_count
xlabel('\theta')
ylabel('$\bar{\psi}$','interpreter','latex')


%%
output_psi_ej_bar=interp1(psi_scale,psi_norm_scale,psi_ej);

% plot(theta_ej,output_psi_ej_bar,'m.');

    
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