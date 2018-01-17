
psiH_diff_PR_map_raw=zeros(NP,Nradial);

for (p=1:NP)
    matrix_size=Nradial;
    
    psi_star_vector=zeros(matrix_size,1);
    psi_star_vector(1)=(psiH_PR_map(p,2)-psiH_PR_map(p,1))/(dist_surf_PR_map(p,2));  
    psi_star_vector(2)=(psiH_PR_map(p,3)-psiH_PR_map(p,1))/(dist_surf_PR_map(p,2)+dist_surf_PR_map(p,3));  
 
    for (r_line=3:matrix_size-1)
        psi_star_vector(r_line)=(psiH_PR_map(p,r_line+1)-psiH_PR_map(p,r_line-1))+(psiH_PR_map(p,r_line)-psiH_PR_map(p,r_line-2));
    end
    psi_star_vector(matrix_size)=(psiH_PR_map(p,matrix_size)-psiH_PR_map(p,matrix_size-1))/(dist_surf_PR_map(p,matrix_size));  
    %psi_star_vector=psi_star_vector';
    
    DR_matrix=zeros(matrix_size,matrix_size);
    %DR_matrix(1,1)=dist_surf_PR_map(2);
    %DR_matrix(1,2)=dist_surf_PR_map(2);
    DR_matrix(1,1)=1;
    DR_matrix(2,2)=1;
    
    frac_diff=0.5;

    for (r_line=3:matrix_size-1)
        DR_matrix(r_line,r_line-2)=frac_diff*dist_surf_PR_map(p,r_line-1);
        DR_matrix(r_line,r_line-1)=(1-frac_diff)*dist_surf_PR_map(p,r_line-1)+dist_surf_PR_map(p,r_line);
        DR_matrix(r_line,r_line)=(dist_surf_PR_map(p,r_line)+(1-frac_diff)*dist_surf_PR_map(p,r_line+1));
        DR_matrix(r_line,r_line+1)=frac_diff*dist_surf_PR_map(p,r_line+1);
    end
    
    DR_matrix(matrix_size,matrix_size)=1;
    
    psiH_diff_PR_map_raw(p,:)=(pinv(DR_matrix)*psi_star_vector)';

  

end


psiH_diff_PR_map=psiH_diff_PR_map_raw;

filt_value=0.95;
filt_value2=1.8;

for (p=1:NP)
    for (r=3:Nradial-3)
        psiH_diff_PR_map(p,r)=((1-filt_value)*psiH_diff_PR_map_raw(p,r-2)+filt_value*psiH_diff_PR_map_raw(p,r-1)+filt_value2*psiH_diff_PR_map_raw(p,r)+filt_value*psiH_diff_PR_map_raw(p,r+1)+(1-filt_value)*psiH_diff_PR_map_raw(p,r+2))/(2+filt_value2);
    end
end
psiH_diff_PR_map_raw=psiH_diff_PR_map;

psiH_diff_PR_map_comp=psiH_diff_PR_map;

for (p=1:NP)
    for (r=2:Nradial-1)
        psiH_diff_PR_map_comp(p,r)=(psiH_PR_map(p,r+1)-psiH_PR_map(p,r-1))/(dist_surf_PR_map(p,r)+dist_surf_PR_map(p,r+1));
    end
end
