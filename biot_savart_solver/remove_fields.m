function struct=remove_fields(struct,fnames,exfnames)
%remove_nonsim_maps Removes maps that are not related to the simulation
%   Detailed explanation goes here
%% Remove unused fields
% Declare which fields need to be removed

if nargin~=3
    exfnames={}; % Excluded names
end

if isempty(fnames)
    fnames=fieldnames(struct);
end

% Complex, but short and effective removal of struct-fields, with check if field .
for j=1:length(fnames)
    if isfield(struct,fnames{j}) && ~any(strcmp(fnames{j},exfnames))
        struct=rmfield(struct,fnames{j});
    end
end

end

