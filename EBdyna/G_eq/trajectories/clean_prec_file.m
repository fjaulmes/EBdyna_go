function [ output_args ] = clean_prec_file( FILENAME,OUTPUTFILENAME )

input=struct();
prec=struct();
load([FILENAME '.mat']);   

par.N_job=input.N_job;

if exist('pop')
prec.pop=pop;
end
if isfield(input,'vpll')
	input.vpll_end=input.vpll;
end
input=remove_fields(input,{'v_end','x_gc','v','Fc_field','x_end','Ekin_end','N_job','vpll'});              % Remove these fields...
if exist('wb')
prec.wb=wb;
prec.wd=wd;
prec.psi_avg=psi_avg;
end
if exist('vpll_avg')
    prec.vpll_avg=vpll_avg;
end
if exist('delta_r')
	prec.delta_r=delta_r;
end
prec.ejected=ejected;


load([OUTPUTFILENAME '.mat'],'output');    
if isfield(output,'x_gc')
	input.x_gc_end=squeeze(output.x_gc(:,:,end));
end
output_prec=remove_fields(output,{'x','v','pphi_kin'});           % Remove these fields...


save([FILENAME '_fix.mat'],'par','input','prec','output_prec');

end

