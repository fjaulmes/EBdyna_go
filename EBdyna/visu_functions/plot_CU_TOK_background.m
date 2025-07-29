% plot a figure with background graphics of this pulse
% PSI_LIMIT=0.0055

OLD_WALL=0;

% first load data
DATA_FOLDER=[SCEN_SIM_FOLDER,'/data_plasma/'];
WALL_FOLDER='./data_common/compassu/';


load([DATA_FOLDER 'XZsmall_fields_tokamak_pre_collapse.mat']);
load([DATA_FOLDER 'motions_map_dimensions.mat']);
try
load([WALL_FOLDER 'wallRZ_CU.mat']);
end

% draw the information
if isempty('gca')
figure
end
hold on
grid on
set(gca,'fontsize',22)
set(gca,'linewidth',4)

ZOFF_FIX=1e22;

[Cvals,ch] =contour(scale_X+R0_grid,scale_Z,psi_XZ',15,'color',[0.4 0 0],'linewidth',1.5);
% contour(scale_X+R0_grid,scale_Z,psi_XZsmall_map',50,'linewidth',1.0);
%# change the ZData property of the inner patches
hh = get(ch,'Children');    %# get handles to patch objects
for i=1:numel(hh)
    zdata = ones(size( get(hh(i),'XData') ));
    set(hh(i), 'ZData',zdata+ZOFF_FIX)
end

PSI_LIMIT=psi_scale(end);
[Cvals,ch]=contour(scale_X+R0_grid,scale_Z,psi_XZ',[PSI_LIMIT PSI_LIMIT],'color',[1.0 0.2 0.2],'linewidth',2.5);
hh = get(ch,'Children');    %# get handles to patch objects
for i=1:numel(hh)
    zdata = ones(size( get(hh(i),'XData') ));
    set(hh(i), 'ZData',zdata+ZOFF_FIX)
end

ph=plot(R_axis,Z_axis,'rx');
ph.ZData=ph.ZData+ZOFF_FIX;

if ~OLD_WALL
    ph=plot(wall_CU.R,wall_CU.Z,'k','linewidth',3);
    ph.ZData=ph.ZData+ZOFF_FIX;
else
    [~ , ~ , ~ , ~ , ph]=plot_CU_vessel(1);
    ph.ZData=ph.ZData+ZOFF_FIX;
    
    hold on
end

axis xy square equal


