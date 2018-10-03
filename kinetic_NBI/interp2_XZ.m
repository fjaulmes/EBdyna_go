function interp_values=interp2_XZ(interp_x,interp_z,FXZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,varargin)
minInputs = 7;
maxInputs = 8;
narginchk(minInputs,maxInputs)

if isempty(varargin)
    PART_LIST=(1:max(size(INDEX_LIST_1)));
else
    PART_LIST=varargin{1};
end
interp_values=(1-interp_x(PART_LIST)).*FXZ_map(INDEX_LIST_1(PART_LIST)).*(1-interp_z(PART_LIST))+(interp_x(PART_LIST)).*FXZ_map(INDEX_LIST_3(PART_LIST)).*(1-interp_z(PART_LIST))+(1-interp_x(PART_LIST)).*FXZ_map(INDEX_LIST_2(PART_LIST)).*(interp_z(PART_LIST))+(interp_x(PART_LIST)).*FXZ_map(INDEX_LIST_4(PART_LIST)).*(interp_z(PART_LIST));
