% plot a figure with background graphics of this pulse
% PSI_LIMIT=0.0055


% first load data
DATA_FOLDER='../../data_tokamak/';
WALL_FOLDER='C:\Users\fabien\EBdyna_go\data_common\';
try
load([DATA_FOLDER 'physics_constants.mat']);
catch
 DATA_FOLDER='./data_tokamak/';
%  WALL_FOLDER='../data_common/';
end
try
load([DATA_FOLDER 'physics_constants.mat']);
catch
 DATA_FOLDER='../data_tokamak/';
%  WALL_FOLDER='../../data_common/'
end
load([DATA_FOLDER 'XZsmall_fields_tokamak_pre_collapse.mat']);
load([DATA_FOLDER 'motions_map_dimensions.mat']);
try
WALL_FOLDER='C:\Users\fabien\EBdyna_go\data_common\';
load([WALL_FOLDER 'MASTU_WALL.mat']);
end

% draw the information
if isempty('gca')
figure
end
hold on
grid on
set(gca,'fontsize',22)
set(gca,'linewidth',4)

contour(scale_X+R0,scale_Z,psi_XZsmall_map',15,'k','linewidth',2);

PSI_LIMIT=psi_scale(end);
contour(scale_X+R0,scale_Z,psi_XZsmall_map',[PSI_LIMIT PSI_LIMIT],'r','linewidth',5);
plot(Raxis,Z_axis,'rx')

plot(MASTU_WALL.R,MASTU_WALL.Z,'k','linewidth',5);

axis xy square equal


