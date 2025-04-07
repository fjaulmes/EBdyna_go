function [struct]=initialize_folder_names_struct(struct)

struct.FINESSE_FOLDER='../';
struct.DATA_FOLDER='../data_tokamak/';
struct.INITIAL_EQUILIBRIUM='../rebuild_equilibrium/';
struct.RECONNECTION_FOLDER='../reconnection_mapping/';
struct.OLLAPSE_MAPS_FOLDER='../build_collapse_maps/';
struct.EPOT_FOLDER='../calculate_Epot/';
struct.EPOT_MAPS_FOLDER='../calculate_Epot/reconnection_maps/';
struct.MOTION_FOLDER='../particles_collapse/';
struct.EQ_FOLDER='../particles_equilibrium/';
end