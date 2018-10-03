function [ leftover ] = lap_psi_fun(P_fac,lap_star_XZ_map_psi,FFprime_XZ_map,Pprime_XZ_map)
leftover =  sum(sum(abs((lap_star_XZ_map_psi-(-FFprime_XZ_map-P_fac*Pprime_XZ_map)))));

end

