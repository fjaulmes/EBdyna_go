function [ TF_coils ,field] = BS_AUG_TFR_apply_coils(TF_map,TF_coils )
%BS_AUG_TFR_apply_coils applies toroidal mode
% Asks for the sign from the add signs function. Then adds superpositions
% the B and A fields of the RMP coils. 

narginchk(2,2);

%% Add them with proper sign
fields_determined=fieldnames(TF_map.TF_coil); % type of fields that have been calculated
for field_type=1:length(fields_determined) % for each field (BX,BY,BZ etc)
    % Pre-allocate
    field.(fields_determined{field_type})=zeros(size(TF_map.TF_coil(1).(fields_determined{field_type})));
    
    for i=1:length(TF_map.TF_coil)
        % Add to field
        field.(fields_determined{field_type})=field.(fields_determined{field_type})...
            +TF_map.TF_coil(i).(fields_determined{field_type});
    end
end

return
end