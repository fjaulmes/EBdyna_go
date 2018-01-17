% *************************************************************
% Mapping the shape of flux surfaces (upper part)
% *************************************************************
%

Z_psi=zeros(Nradial,NR)+Z_axis;
Z_psi_fit=zeros(Nradial,NR)+Z_axis;





% finding what this the inner (left) and outer (right)
% boundaries for the external flux surface
[eps x_inner ] = max(psi2D_mask(mid_X:-1:1,Z_axis_pos));
x_inner=mid_X-x_inner;
[eps x_outer ] = max(psi2D_mask(mid_X:NR,Z_axis_pos));
x_outer=x_outer+mid_X;

Z_psi(Nradial,x_inner)=0;
Z_psi(Nradial,x_outer)=0;





disp('... Building upper part of Z_psi curves ...');

% the first two core psi values are ignored
FITSTEP=1;
%FITSTEP=2+ASDEX_LIKE_EQUILIBRIUM;
%for (n=2:4)
    %fit_flux_up_fitstep;
	%Z_psi(n,X1_Nradial(n):X2_Nradial(n))=Z_contour_up(n,1:length_contour_up(n));
%end
%Z_psi(2,X1_Nradial(2):X2_Nradial(2))=Z_contour_up(2,1:length_contour_up(2));

err_value=1;
for (n=2:Nradial-4)
	try
		fit_flux_up_fitstep;
		err_value=0; 
	catch err
        err_value=1;
	end
	while err_value==1
        err_value=0; 	    
		FITSTEP=FITSTEP+1;
		n
		try
			fit_flux_up_fitstep;
		catch err
			err_value=1;
	    end
	end
end


err_value=1;
for (n=Nradial-3:Nradial-1)
	try
		Z_psi(n,X1_Nradial(n):X2_Nradial(n))=interp1(X_contour_up(n,FITSTEP:FITSTEP:length_contour_up(n)-FITSTEP),Z_contour_up(n,FITSTEP:FITSTEP:length_contour_up(n)-FITSTEP),X_scale(X1_Nradial(n):X2_Nradial(n)));
		err_value=0; 
	catch err
        err_value=1;
	end
	while (err_value==1) && (FITSTEP<Nradial)		
        err_value=0;		
		FITSTEP=FITSTEP+1;
		n
		try
			Z_psi(n,X1_Nradial(n):X2_Nradial(n))=interp1(X_contour_up(n,FITSTEP:FITSTEP:length_contour_up(n)-FITSTEP),Z_contour_up(n,FITSTEP:FITSTEP:length_contour_up(n)-FITSTEP),X_scale(X1_Nradial(n):X2_Nradial(n)));
		catch err
			err_value=1;
	    end
	end
end
%if FITSTEP>(Nradial/20)
Z_psi(Nradial,:)=Z_psi(Nradial-1,:);
%end




%A(isnan(A))=0

Z_psi(isnan(Z_psi))=Z_axis;

%Nradial_min=6;

Nradial_inf=min(find(X2_Nradial-X1_Nradial>2));

for (n=3:Nradial)
    Z_psi_fit(n,X1_Nradial(n):X2_Nradial(n))=(interp1(DX*[X1_Nradial(n):X2_Nradial(n)],Z_psi(n,X1_Nradial(n):X2_Nradial(n)),DX*[X1_Nradial(n):X2_Nradial(n)],'cubic'));
    Z_psi(n,:)=(Z_psi_fit(n,:)+Z_psi(n,:))/2;
    %Z_psi(n,X1(n):X2(n))=(interp1(DX*[X1(n):8:X2(n)],Z_psi(n,X1(n):8:X2(n)),DX*[X1(n):X2(n)],'cubic'));
end
    
Z_psi=max(Z_psi,Z_axis);


for (n=2:Nradial)
    Z_psi_fit(n,X1_Nradial(n):X2_Nradial(n))=Z_psi(n,X1_Nradial(n):X2_Nradial(n));  
end

Z_psi_fit(isnan(Z_psi_fit))=Z_axis;

Z_psi_fit(1,:)=zeros(1,NR)+Z_axis;





for (n=2:Nradial)
    Z_psi_fit(n,1:X1_Nradial(n))=Z_axis;
    Z_psi_fit(n,X2_Nradial(n):NR)=Z_axis;

    
    Z_psi(n,1:X1_Nradial(n))=Z_axis;
    Z_psi(n,X2_Nradial(n):NR)=Z_axis;

end
Z_psi_fit=0.5*(Z_psi_fit+Z_psi);
Z_psi_fit=max(Z_psi_fit,Z_axis);

Z_psi_fit_up=Z_psi_fit;

% Display the result of the polynomial interpolation

figure(2);
hold on

axis xy
xlabel('R (m)');
ylabel('Z (m)');
[ny, nx] = size(psi2D'); 

for (n=1:2:Nradial)
    plot(X_scale,Z_psi_fit(n,:),'g');
%     plot(xpos(X1_Nradial(n):X2_Nradial(n)),Z_psi_fit(n,X1_Nradial(n):X2_Nradial(n)),'k');
end

% xlim([Rpos(X1(n)) Rpos(X2(n))]);
% ylim([zpos(mid_Z) zpos(sup_Z)]);




%disp('... calculate_ki_angle_geometry.m ...');

%run('calculate_ki_angle_geometry');

