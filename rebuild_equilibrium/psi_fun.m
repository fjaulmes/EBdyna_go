function [ leftover ] = psi_fun(psi_fac,BZ,BZ_XZ_map)
leftover =  sum(sum(abs((psi_fac*BZ-BZ_XZ_map))));

end

