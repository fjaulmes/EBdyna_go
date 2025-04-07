% Markers selection
load('./output/G_eq_1064_full_60k.mat');

input_new=input;
input_new.v=input_new.v_end;
input_new.x=input_new.x_end;
input_new=remove_fields(input_new,{'vpll_ini','mm','pphi_kin','x_ini','v_ini','x_end','v_end','x_gc_end','N_job','particle_nr'});
fnames=fieldnames(input_new);


%     vpll_ini: [49759x1 double]
%           mm: [49759x1 double]
%     pphi_kin: [49759x1 double]
    
for j=1:length(fnames)
    fname=fnames{j};
    data_input=input_new.(fname);
    if size(data_input,1)>1
        if size(data_input,2)==1
    input_new.(fname)=data_input(~ejected);
        end
        input_new.(fname)=data_input(~ejected,:);
    else
        input_new.(fname)=data_input;
    end
end

input_new.N_total=length(input_new.Ekin);
input=input_new

% Here we sace the pure sd distribution
% for prec statistics calculations
save('./input/G_eq_1064_sd15_50k.mat','input')


disp('OFFSET FOR BIRTH MATRIX TO ENTER IN SIM PARAMETERS:');
disp(input.N_total);
input.SD_MARKERS_END=input.N_total;

% recycling of the original distribution of NBI ionization markers
% this is a bit ugly statistics but saves us much trouble to re-initialize
% markers when doing a simulation with sd and NBI source
% note that this can be replace by any list of relevant ionization markers

% old_dist=load('../input/BBNBI_R0p65_incl5_60k.mat')

par.LOADNAME='./input/BBNBI_R0p65_incl5_60k.mat';


% input_ionization=old_dist.input;

par.NB_PROCESS=1;
par.PROCESS_NUMBER=1;
[maps dim]=initialize_maps(1,0);
[x,v,input_ionization,output_ionization,ejected_ionization]=initialize_particles();

input.x=[input.x ; x];
input.v=[input.v ; v];
input.Ekin=[input.Ekin ; input_ionization.Ekin];
input.N_total=length(input.Ekin);

save('./input/G_eq_1064_sd15_50k60k.mat','input')
input
