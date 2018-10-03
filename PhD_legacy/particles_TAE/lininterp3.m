function v = lininterp3( V,INDEX_LIST3D_1,INDEX_LIST3D_2,INDEX_LIST3D_3,INDEX_LIST3D_4,INDEX_LIST3D_5,INDEX_LIST3D_6,INDEX_LIST3D_7,INDEX_LIST3D_8,slopex,slopey,slopez)
% linear interpolation, given set of X, Y, Z, and V values, and an x, y, z query
% assumes X, Y, and Z values are in strictly increasing order
%
% Differences from matlab built-in :
%       order of arguments switched
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%

%if ((length(X) ~= size(V, 1)) || (length(Y) ~= size(V, 2)) || (length(Z) ~= size(V, 3))), error('[length(X), length(Y), length(Z)] does not match size(V)'); end

% v=zeros(length(x),1);
% [INDEX_LIST3D_1 INDEX_LIST3D_2 INDEX_LIST3D_3 INDEX_LIST3D_4 INDEX_LIST3D_5 INDEX_LIST3D_6 INDEX_LIST3D_7 INDEX_LIST3D_8 slopex slopey slopez] = build_3Dinterp_indexarrays(X, Y, Z, x, y, z);



    v = V(INDEX_LIST3D_1) .* (1 - slopex) .* (1 - slopey) .* (1 - slopez) + V(INDEX_LIST3D_2) .* slopex .* (1 - slopey) .* (1 - slopez) ...
      + V(INDEX_LIST3D_3) .* (1 - slopex) .*       slopey .* (1 - slopez) + V(INDEX_LIST3D_4) .* slopex .*       slopey .* (1 - slopez) ...
      + V(INDEX_LIST3D_5) .* (1 - slopex) .* (1 - slopey) .*       slopez + V(INDEX_LIST3D_6) .* slopex .* (1 - slopey) .*       slopez ...
      + V(INDEX_LIST3D_7) .* (1 - slopex) .*       slopey .*       slopez + V(INDEX_LIST3D_8) .* slopex .*       slopey .*       slopez ;

end
