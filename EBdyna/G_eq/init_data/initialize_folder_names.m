function folders=intialize_folder_names
global par


% the given folders are relative to the G_eq/.. (EBdyna) folder

DATA_SHOT=[par.sim_folder_string,'/data_plasma/'];
DATA_SHOT_REF=['./data_plasma/',par.tokamak,'/',par.shot_name,'/'];
GEQ_INPUT=[par.sim_folder_string,'/input/'];
GEQ_OUTPUT=[par.sim_folder_string,'/output/'];
DATA_COMMON_TOKAMAK=['./data_common/',par.tokamak,'/'];
DATA_COMMON_PHYSICS=['./data_common/physics_data/'];

if ~isfolder(GEQ_OUTPUT)
    mkdir(GEQ_OUTPUT);
end


folders=struct();
folders.DATA_SHOT=DATA_SHOT;
folders.DATA_SHOT_REF=DATA_SHOT_REF;
folders.GEQ_INPUT=GEQ_INPUT;
folders.GEQ_OUTPUT=GEQ_OUTPUT;
folders.DATA_COMMON_TOKAMAK=DATA_COMMON_TOKAMAK;
folders.DATA_COMMON_PHYSICS=DATA_COMMON_PHYSICS;


par.folders=folders;
