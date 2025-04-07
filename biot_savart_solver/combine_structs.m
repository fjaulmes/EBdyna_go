function [ struct1 ] = combine_structs( struct1,struct2)
%combine_structs Combines 2 structs if no fields with the same name are
%present
%   Is usefull for individual load commands. Loops over the fieldnames of
%   struct 2 and adds them to struct 1.
fnames=fieldnames(struct2);
for i=1:length(fnames)
    if isfield(struct1,fnames{i})
        if isequal(struct1.(fnames{i}),struct2.(fnames{i}))
            warning('Combining 2 structs with duplicate data')
            continue
        end
        error('Field with same name in both structs')
    end
    struct1.(fnames{i})=struct2.(fnames{i});
end

end