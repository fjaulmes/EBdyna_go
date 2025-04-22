function [n0,n0D2] = neutrals_2D_interp(x)
    %NEUTRALS_2D_INTERP query the values of neutral density in 2D 
    %   For more accurate CX calculations
    global par maps dim
    
    
    % Indexes 2D
    X_ind=((x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
    Z_ind=( x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
    
    %% BA_INTERP2 LINEAR
    n0=ba_interp2(maps(1).neutrals_XZ,Z_ind,X_ind,'linear');
    if par.N0_FAC_D2>0
        n0D2=ba_interp2(maps(1).neutralsD2_XZ,Z_ind,X_ind,'linear');
    else
        n0D2=[];
    end
end

