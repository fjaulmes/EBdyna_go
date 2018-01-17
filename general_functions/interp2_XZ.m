function values=interp2_XZ(interp_x,interp_z,map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,expr)
% minInputs = 7;
% maxInputs = 8;
% narginchk(minInputs,maxInputs)

%% Nargin check
if nargin <3 || nargin==5 || nargin==6 || nargin>8
    error('Incorrect number of input arguments')
elseif nargin==3 || nargin ==4
    %% New structured EBdyna_go
    % Rename to proper variables
    indexes = interp_x;
    slopes = interp_z;
    if nargin==3
        expr=true(numel(slopes.x),1);
    elseif nargin==4
        expr=INDEX_LIST_1;
    end
    % Linear interpolation
    sxsz=slopes.x(expr).*slopes.z(expr);
    values=...
         (1-slopes.x(expr)-slopes.z(expr)+sxsz  ).*map(indexes.lb(expr))...   %left  bottom corner
        +(                 slopes.z(expr)-sxsz  ).*map(indexes.lt(expr))...   %left  top    corner
        +(  slopes.x(expr)-sxsz                 ).*map(indexes.rb(expr))...   %right bottom corner
        +(   sxsz                               ).*map(indexes.rt(expr));     %right top    corner
    return
end

%% Support for old EBdyna_go
if nargin==7
    expr=true(size(INDEX_LIST_1));
end
values=(1-interp_x(expr)).*map(INDEX_LIST_1(expr)).*(1-interp_z(expr))+(interp_x(expr)).*map(INDEX_LIST_3(expr)).*(1-interp_z(expr))+(1-interp_x(expr)).*map(INDEX_LIST_2(expr)).*(interp_z(expr))+(interp_x(expr)).*map(INDEX_LIST_4(expr)).*(interp_z(expr));
end
