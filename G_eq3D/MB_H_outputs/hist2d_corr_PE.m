function [ data_hist ] = Untitled2(pos_R_gc, pos_Z_gc, R_BINS, Z_BINS, part_weight,part_weight_Ekin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
histR_size=length(R_BINS)-1;
histZ_size=length(Z_BINS)-1;

data_hist=zeros(histR_size,histZ_size);

for Rb=1:histR_size
    for Zb=1:histZ_size
        PART_POL_SECTION=find((pos_R_gc>=R_BINS(Rb)).*(pos_R_gc<R_BINS(Rb+1)).*(pos_Z_gc>=Z_BINS(Zb)).*(pos_Z_gc<Z_BINS(Zb+1)));
        data_hist(Rb,Zb)=sum(part_weight(PART_POL_SECTION).*part_weight_Ekin(PART_POL_SECTION));
    end
end


end

