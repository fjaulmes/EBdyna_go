%% code to calculate Janev based cross section

function sigma=sigma_Janev(E,ne,Te,Zeff,A,B)
    
lnN=log(ne/1e19);
E=log(E);
U=log(Te);


for i=1:3
    for j=1:3
        for k=1:3
        Em(i,j,k)=Em^(i-1);
        lnNm(i,j,k)=lnN^(j-1);
        U(i,j,k)=U^(k-1);
        end
    end
end

sigH=1e-16/E*exp(sum(A.*Em.*lnNm.*U,'all'));
Sz=sum(B.*Em.*lnNm.*U,'all');
sigma=sigH(1+(Zeff-1)*Sz);

end



