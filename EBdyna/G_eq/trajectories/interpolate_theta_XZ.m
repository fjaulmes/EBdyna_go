function [theta] = interpolate_theta_XZ(indexes,slopes,expr_calc)
global maps par
persistent T_N T_PS LAST_SCHEME

if isfield(maps,'theta') && isstruct(maps(1).theta) 
    %% Older algorithm (where the two theta maps are in a struct)
	% Expr_calc
	if nargin==2
		expr_calc=true(size(indexes.lb));
		theta=maps(1).theta.normal_XZ(indexes.lb); % Estimate of theta to decide map
	else
		theta=NaN(size(expr_calc));
		theta(expr_calc)=maps(1).theta.normal_XZ(indexes.lb(expr_calc)); % Estimate of theta to decide map
	end
	
	% Find expression for which map
	expr_normal=theta>0.5*pi & theta<1.5*pi & expr_calc; % Normal map 
	expr_phase_shift=~expr_normal           & expr_calc;              % Map with phase difference
	
	theta(expr_normal     )=interp2_XZ(indexes,slopes,maps(1).theta.normal_XZ     ,expr_normal);
	theta(expr_phase_shift)=interp2_XZ(indexes,slopes,maps(1).theta.phase_shift_XZ,expr_phase_shift);
	
	% Make sure theta of phase shifted are in domain [0,2*pi]
	theta(expr_phase_shift)=mod(theta(expr_phase_shift),2*pi);
    
else
    %% New algorithm - faster and theta's not a substruct
    if par.interp_scheme==1
        indexes.lb(isnan(indexes.lb))=1;
        theta=maps(1).theta_normal_XZ(indexes.lb);
        if nargin==2
            expr_calc=true(size(theta));
        end
    else
         ind_X=indexes;
         ind_Z=slopes;
         theta=ba_interp2(maps(1).theta_normal_XZ,ind_Z,ind_X,'nearest');
    end
    
	% Find expression for which map
	expr_normal=theta>0.5*pi & theta<1.5*pi; % Normal map 
	expr_phase_shift=~expr_normal          ;              % Map with phase difference
	
    switch par.interp_scheme
        case 1
            theta(expr_normal      & expr_calc)=interp2_XZ(indexes,slopes,maps(1).theta_normal_XZ     ,expr_normal & expr_calc);
            theta(expr_phase_shift & expr_calc)=interp2_XZ(indexes,slopes,maps(1).theta_phase_shift_XZ,expr_phase_shift & expr_calc);
        case 2
            theta(expr_normal)     =interp2(maps(1).theta_normal_XZ     ,ind_Z(expr_normal)     ,ind_X(expr_normal)     ,'*linear');
            theta(expr_phase_shift)=interp2(maps(1).theta_phase_shift_XZ,ind_Z(expr_phase_shift),ind_X(expr_phase_shift),'*linear');
        case 3
            theta(expr_normal)     =ba_interp2(maps(1).theta_normal_XZ     ,ind_Z(expr_normal)     ,ind_X(expr_normal)     ,'linear');
            theta(expr_phase_shift)=ba_interp2(maps(1).theta_phase_shift_XZ,ind_Z(expr_phase_shift),ind_X(expr_phase_shift),'linear');
        case 4
            if isempty(T_N) || LAST_SCHEME~=par.interp_scheme
                T_N =griddedInterpolant(maps(1).theta_normal_XZ,'linear');
                T_PS=griddedInterpolant(maps(1).theta_phase_shift_XZ,'linear');
            end
            theta(expr_normal)      =T_N (ind_X(expr_normal)        ,ind_Z(expr_normal));
            theta(expr_phase_shift) =T_PS(ind_X(expr_phase_shift)   ,ind_Z(expr_phase_shift));
        case 5
            theta(expr_normal)     =interp2(maps(1).theta_normal_XZ     ,ind_Z(expr_normal)     ,ind_X(expr_normal)     ,'*cubic');
            theta(expr_phase_shift)=interp2(maps(1).theta_phase_shift_XZ,ind_Z(expr_phase_shift),ind_X(expr_phase_shift),'*cubic');
        case 6
            theta(expr_normal)     =ba_interp2(maps(1).theta_normal_XZ     ,ind_Z(expr_normal)     ,ind_X(expr_normal)     ,'cubic');
            theta(expr_phase_shift)=ba_interp2(maps(1).theta_phase_shift_XZ,ind_Z(expr_phase_shift),ind_X(expr_phase_shift),'cubic');
        case 7
            if isempty(T_N) || LAST_SCHEME~=par.interp_scheme
                T_N=griddedInterpolant(maps(1).theta_normal_XZ,'cubic');
                T_PS=griddedInterpolant(maps(1).theta_phase_shift_XZ,'cubic');
            end
            theta(expr_normal)      =T_N (ind_X(expr_normal)        ,ind_Z(expr_normal));
            theta(expr_phase_shift) =T_PS(ind_X(expr_phase_shift)   ,ind_Z(expr_phase_shift));
        case 8
            theta(expr_normal)     =interp2(maps(1).theta_normal_XZ     ,ind_Z(expr_normal)     ,ind_X(expr_normal)     ,'*spline');
            theta(expr_phase_shift)=interp2(maps(1).theta_phase_shift_XZ,ind_Z(expr_phase_shift),ind_X(expr_phase_shift),'*spline');
        case 9
            if isempty(T_N) || LAST_SCHEME~=par.interp_scheme
                T_N=griddedInterpolant(maps(1).theta_normal_XZ,'spline');
                T_PS=griddedInterpolant(maps(1).theta_phase_shift_XZ,'spline');
            end
            theta(expr_normal)      =T_N (ind_X(expr_normal)        ,ind_Z(expr_normal));
            theta(expr_phase_shift) =T_PS(ind_X(expr_phase_shift)   ,ind_Z(expr_phase_shift));
    end
    % Make sure theta of phase shifted are in domain [0,2*pi]
    theta(theta<0)=theta(theta<0)+2*pi;         % Faster then mod-function
    LAST_SCHEME=par.interp_scheme;
end