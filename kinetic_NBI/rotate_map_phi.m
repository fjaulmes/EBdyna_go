function [ output_map ] = rotate_map_phi( input_map,phi_rank )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        size_radial=size(input_map,1);
        size_omega=size(input_map,2);
        omega_phi_rank=max(size_omega-(phi_rank),1);
        output_map(:,:)=[input_map(:,omega_phi_rank:size_omega)  input_map(:,2:omega_phi_rank)];
%         output_map=output_map';

end

