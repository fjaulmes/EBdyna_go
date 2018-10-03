% *************************************************************
% Discretisation to a finite number of flux values
% *************************************************************
%
disp('fitting_Zpsi_down_profiles.m ...');


%NPSI=round(NX*0.18)+1;
Z_psi_down=zeros(Nradial,NR)+Z_axis;
Z_psi_fit_down=zeros(Nradial,NR)+Z_axis;

%n=2;
%FITSTEP=2;
%Z_psi_down(n,X1_Nradial(n):X2_Nradial(n))=interp1(X_contour_down(n,2:FITSTEP:length_contour_down(n)),Z_contour_down(n,2:FITSTEP:length_contour_down(n)),X_scale(X1_Nradial(n):X2_Nradial(n)));

%FITSTEP=1+ASDEX_LIKE_EQUILIBRIUM;
disp('... Building lower part of Z_psi curves  ...');

%for (n=2:5)
%    fit_flux_down_fitstep;
	%Z_psi(n,X1_Nradial(n):X2_Nradial(n))=Z_contour_down(n,length_contour_down(n));
%end
Z_psi(2,X1_Nradial(2):X2_Nradial(2))=Z_contour_down(2,length_contour_down(2));

FITSTEP=1;

err_value=1;
for (n=6:Nradial-4)
	try
		fit_flux_down_fitstep;
		err_value=0; 
	catch err
        err_value=1;
	end
	while err_value==1
        err_value=0; 	    
		FITSTEP=FITSTEP+1;
		n
		try
			fit_flux_down_fitstep;
		catch err
			err_value=1;
	    end
	end
end

err_value=1;
for (n=Nradial-3:Nradial-1)
	try
		Z_psi_down(n,X1_Nradial(n):X2_Nradial(n))=interp1(X_contour_down(n,FITSTEP:FITSTEP:length_contour_down(n)-FITSTEP),Z_contour_down(n,FITSTEP:FITSTEP:length_contour_down(n)-FITSTEP),X_scale(X1_Nradial(n):X2_Nradial(n)));
		err_value=0; 
	catch err
        err_value=1;
	end
	while (err_value==1) && (FITSTEP<Nradial)		
        err_value=0;		
		FITSTEP=FITSTEP+1;
		n
		try
			Z_psi_down(n,X1_Nradial(n):X2_Nradial(n))=interp1(X_contour_down(n,FITSTEP:FITSTEP:length_contour_down(n)-FITSTEP),Z_contour_down(n,FITSTEP:FITSTEP:length_contour_down(n)-FITSTEP),X_scale(X1_Nradial(n):X2_Nradial(n)));
		catch err
			err_value=1;
	    end
	end
end
%if FITSTEP>(Nradial/20)
Z_psi_down(Nradial,:)=Z_psi_down(Nradial-1,:);
%end


Z_psi_down=min(Z_psi_down,Z_axis);

NP_half=round(NP/2)+1;

%Nradial_min=6;
% 
Nradial_inf=min(find(X2_Nradial-X1_Nradial>2));

for (n=Nradial_inf:Nradial-1)
    Z_psi_fit_down(n,X1_Nradial(n):X2_Nradial(n))=(interp1(DX*[X1_Nradial(n):2:X2_Nradial(n)],Z_psi_down(n,X1_Nradial(n):2:X2_Nradial(n)),DX*[X1_Nradial(n):X2_Nradial(n)],'cubic'));
    Z_psi_down(n,:)=(Z_psi_fit_down(n,:)+Z_psi_down(n,:))/2;
end

Z_psi_down(isnan(Z_psi_down))=Z_axis;
Z_psi_down=min(Z_psi_down,Z_axis);

% for (n=Nradial_min:Nradial)
%     %disp('n value = ');disp(n);
%     Raxis_shift=R0-R_inf_lim+xi_psi_Nradial(n);
%     n_100=(n/Nradial)*100;
%     pol_deg=round((n_100)/20+2);
%     Z_psi_pol_1=polyfit(DX*([X1_Nradial(n):axis_pos(n)]-1)-Raxis_shift,Z_psi_down(n,X1_Nradial(n):axis_pos(n)),pol_deg);
%     Z_psi_pol_top=polyfit(DX*([axis_pos(n)-D_flat_axis_Nradial(n):axis_pos(n)+D_flat_axis_Nradial(n)]-1)-Raxis_shift,Z_psi_down(n,axis_pos(n)-D_flat_axis_Nradial(n):axis_pos(n)+D_flat_axis_Nradial(n)),4);
%     Z_psi_pol_2=polyfit(DX*([axis_pos(n):X2_Nradial(n)]-1)-Raxis_shift,Z_psi_down(n,axis_pos(n):X2_Nradial(n)),pol_deg);
%     
%     Z_psi_fit_down(n,X1_Nradial(n):axis_pos(n)-D_flat_axis_Nradial(n))=polyval(Z_psi_pol_1,DX*([X1_Nradial(n):axis_pos(n)-D_flat_axis_Nradial(n)]-1)-Raxis_shift);
%     Z_psi_fit_down(n,axis_pos(n)+D_flat_axis_Nradial(n):X2_Nradial(n))=polyval(Z_psi_pol_2,DX*([axis_pos(n)+D_flat_axis_Nradial(n):X2_Nradial(n)]-1)-Raxis_shift);
%     Z_psi_fit_down(n,axis_pos(n)-D_flat_axis_Nradial(n):axis_pos(n)+D_flat_axis_Nradial(n))=polyval(Z_psi_pol_top,DX*([axis_pos(n)-D_flat_axis_Nradial(n):axis_pos(n)+D_flat_axis_Nradial(n)]-1)-Raxis_shift);
%     Z_psi_fit_down(n,X1_Nradial(n)+1)=Z_psi_down(n,X1_Nradial(n)+1);
%     Z_psi_fit_down(n,X2_Nradial(n)-1)=Z_psi_down(n,X2_Nradial(n)-1);
%     
%     pol_deg=round(0.12*n_100+3);
%     
%     Z_psi_pol=polyfit(2*(DX*([X1_Nradial(n)+round(0.1*n_100):X2_Nradial(n)-round(0.1*n_100)]-1)-Raxis_shift),Z_psi_fit_down(n,X1_Nradial(n)+round(0.1*n_100):X2_Nradial(n)-round(0.1*n_100)),pol_deg);
%     Z_psi_fit_down(n,X1_Nradial(n)+2:X2_Nradial(n)-2)=polyval(Z_psi_pol,2*(DX*([X1_Nradial(n)+2:X2_Nradial(n)-2]-1)-Raxis_shift));
%     Z_psi_fit_down(n,:)=(Z_psi_fit_down(n,:)+Z_psi_down(n,:))/2;
%     Z_psi_fit_down(n,X1_Nradial(n))=0;
%     Z_psi_fit_down(n,X2_Nradial(n))=0;
%     
%     Z_psi_pol=polyfit(2*(DX*([X1_Nradial(n)+1:X2_Nradial(n)-1]-1)-Raxis_shift),Z_psi_fit_down(n,X1_Nradial(n)+1:X2_Nradial(n)-1),pol_deg);
%     Z_psi_fit_down(n,X1_Nradial(n)-2:X2_Nradial(n)+2)=polyval(Z_psi_pol,2*(DX*([X1_Nradial(n)-2:X2_Nradial(n)+2]-1)-Raxis_shift));
%     
% end


Z_psi_fit_down(1,:)=zeros(1,NR)+Z_axis;

for (n=2:Nradial)
    Z_psi_fit_down(n,X1_Nradial(n):X2_Nradial(n))=Z_psi_down(n,X1_Nradial(n):X2_Nradial(n));
end

Z_psi_fit_down(isnan(Z_psi_fit_down))=Z_axis;


for (n=1:Nradial)
    Z_psi_fit_down(n,1:X1_Nradial(n))=Z_axis;
    Z_psi_fit_down(n,X2_Nradial(n):NR)=Z_axis;

    
    Z_psi_down(n,1:X1_Nradial(n))=Z_axis;
    Z_psi_down(n,X2_Nradial(n):NR)=Z_axis;

end

% Z_psi_fit_down=0.5*(Z_psi_fit_down+Z_psi_down);
Z_psi_fit_down=min(Z_psi_fit_down,Z_axis);


% Display the result of the polynomial interpolation

figure(2);
hold on

for (n=1:2:Nradial)
    plot(X_scale,Z_psi_fit_down(n,:),'g');
%     plot(X_scale(X1_Nradial(n):X2_Nradial(n)),Z_psi_fit_down(n,X1_Nradial(n):X2_Nradial(n)),'r');
end


