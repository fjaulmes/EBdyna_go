

SAVE_INI=1
SAVE_FINAL=1


FILENAME_INI='dist_NBI_PRcorr_AUG31557_2p25_EpitchRZ_ini.dat'
FILENAME_END='dist_NBI_PRcorr_AUG31557_2p25_EpitchRZ_125_end.dat'
% 
string_title='# initial NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
% 

% save dist_EpitchRZ_ini.dat -append  string_title -ASCII
if SAVE_INI==1
    save(FILENAME_INI,'-append','Raxis','Z_axis','-ASCII')
    save(FILENAME_INI,'-append','l1','l2','l3','l4','-ASCII')
    save(FILENAME_INI,'-append','Ekin_values','-ASCII')
    save(FILENAME_INI,'-append','pitch_values','-ASCII')
    save(FILENAME_INI,'-append','R_values','-ASCII')
    save(FILENAME_INI,'-append','Z_values','-ASCII')
end

for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        data_array=squeeze(dist_EpitchRZ_ini(eb,pb,:,:));
        data_array=data_array';
        if SAVE_INI==1
            save(FILENAME_INI,'-append','data_array','-ASCII')
        end
    end
end






string_title='# final NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
%save dist_EpitchRZ_end.dat -append string_title -ASCII
if SAVE_FINAL==1
    save(FILENAME_END,'-append','Raxis','Z_axis','-ASCII')
    save(FILENAME_END,'-append','l1','l2','l3','l4','-ASCII')
    save(FILENAME_END,'-append','Ekin_values','-ASCII')
    save(FILENAME_END,'-append','pitch_values','-ASCII')
    save(FILENAME_END,'-append','R_values','-ASCII')
    save(FILENAME_END,'-append','Z_values','-ASCII')
end

for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        data_array=squeeze(dist_EpitchRZ_end(eb,pb,:,:));
        data_array=data_array';

        if SAVE_FINAL==1
            save(FILENAME_END,'-append','data_array','-ASCII')
        end
    end
end


%%

figure(1)
set(gca,'fontsize',20)
contourf(R_values,Z_values,1e6*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(dist_EpitchRZ_end(1:end,:,:,:),2),1))'),12);axis xy;hold on

%%
figure(2)
set(gca,'fontsize',20)
contourf(R_values,Z_values,1e6*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(poloidal_hist_ini(1:end,:,:,:),2),1))'),12);axis xy;hold on


