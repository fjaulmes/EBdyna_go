% describe_contour allows to concatenate pieces of contours that are for the same value
% returns nb_split=-1 if the data is not given in a way that the contour is continuous
% DLMAX gives the tolerance for two points to make a continuous curve
% one can take typically DLMAX=4
function [psi_value Nlength_global nb_pieces contour_data_x contour_data_z] = describe_contour(Nbegin, ContPsiStar, DLMAX)
    MIN_CONTOUR_SIZE=80;
    MAX_NRANK=size(ContPsiStar,2);
	psi_value=ContPsiStar(1,Nbegin);
	psi_value_partial=psi_value;
	Nlength_global=ContPsiStar(2,Nbegin);
	Nlength_partial=Nlength_global;
	n_rank=Nbegin;
	clear contour_data_x;
	clear contour_data_z;
    contour_data_x=zeros(1,1);
    contour_data_z=zeros(1,1);
    contour_data_x(:,:)=[];
    contour_data_z(:,:)=[];
    nb_pieces=0;

	while(psi_value_partial==psi_value)
	    nb_pieces=nb_pieces+1;
		Nlength_partial=ContPsiStar(2,n_rank);
        if (Nlength_partial>MIN_CONTOUR_SIZE)
            if (~isempty(contour_data_x))
                contour_data_x=[contour_data_x ContPsiStar(1,n_rank+1:n_rank+Nlength_partial)];
                contour_data_z=[contour_data_z ContPsiStar(2,n_rank+1:n_rank+Nlength_partial)];
            else
                contour_data_x=[ ContPsiStar(1,n_rank+1:n_rank+Nlength_partial)];
                contour_data_z=[ ContPsiStar(2,n_rank+1:n_rank+Nlength_partial)];
            end
        end
		x_value=ContPsiStar(1,n_rank+Nlength_partial);
		z_value=ContPsiStar(2,n_rank+Nlength_partial);
		if ((n_rank+Nlength_partial)<=MAX_NRANK)
			n_rank=n_rank+Nlength_partial;
			psi_value_partial=ContPsiStar(1,n_rank);
			if (psi_value_partial==psi_value) && (Nlength_partial>MIN_CONTOUR_SIZE)
				Nlength_global=Nlength_global+Nlength_partial;
				if (abs(x_value-ContPsiStar(1,n_rank+2))<DLMAX) && abs(z_value-ContPsiStar(2,n_rank+2))<DLMAX
					nb_pieces=nb_pieces+1;
                else
                    nb_pieces=-1;				
				end
			end
		end
	end