Ekin_vD_indiv_evol=zeros(NB_TIME_STAMPS,1);

for t=2:NB_TIME_STAMPS
    Ekin_vD_indiv_evol(t)=Ekin_vD_indiv_evol(t-1)+part2wave_vD_power_output(t,601)*h;
    Ekin_indiv_evol(t)=Ekin_vD_indiv_evol(t-1)+part2wave_power_output(t,601)*h;
end