clear final_psi_values

final_psi_values(1)=Psih_final_recalc(3);

if (f<=11)||((f>=EXPULSION_FRAME_NUMBER-1)&&(f<=EXPULSION_FRAME_NUMBER+1))
    %higher number of contours required for the expulsion time and init
    final_psi_values=[final_psi_values max(Psih_final_recalc(3:2:Npsi_values),0)];
elseif (f<1001)
    final_psi_values=[final_psi_values max(Psih_final_recalc(5:2:Npsi_values),0)];
else
    final_psi_values=[final_psi_values max(Psih_final_recalc(10:10:Npsi_values),0)];
end
final_psi_values(1)=psi_star_final(1);

% adjust lowest psi values
final_psi_values(end-1)=0.6*final_psi_values(end-1)+0.4*final_psi_values(end-2);
final_psi_values(end)=0.5*final_psi_values(end)+0.5*final_psi_values(end-1);
final_psi_values(end+1)=0.5*final_psi_values(end);

Npsi_values=length(final_psi_values);




% First getting the contours in region 3
% Trying to get as close to separatrix as possible

if (psi_limit13>final_psi_values(end))
    psi_sep_rank=round(interp1(final_psi_values,(1:Npsi_values),psi_limit13));
else
    psi_sep_rank=length(final_psi_values);
end
if (psi_limit13>=abs(Psih_final_recalc(1))-1e-6)
    psi_sep_rank=1
end

Npsi_values_region3=psi_sep_rank;