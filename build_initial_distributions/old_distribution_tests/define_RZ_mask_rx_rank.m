
mask_XZ_zoom_map_reconnection=zeros(NZ_zoom,NZ_zoom);

for (x=minX:maxX)
    for (z=minZ:maxZ)
        if(radial_XZ_zoom_map(x,z)<(rx_rank-delta_rx))
            mask_XZ_zoom_map_reconnection(x,z)=1;
        end
    end
end

% mask_XZ_zoom_reconnection=zeros(4*NX,4*NX);
% 
% for (x=2:4*NX-1)
%     for (z=2:4*NX-1)
%         if(radial_XZ_zoom_map(x,z)<rx_rank-delta_rx)&&(radial_XZ_zoom_map(x,z)~=0)
%             mask_XZ_zoom_reconnection(x,z)=1;
%         end
%     end
% end

