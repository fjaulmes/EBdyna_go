%% my_rand
% build in rand.m funciton in matlab is not inclusive on the interval 0 to
% 1, this is an attempt to correct this. This function will create a
% uniform distribution of peudo-random numbers [0,1[.
%
% See also rand, randi

function r = my_rand(arr)
    p = 1e7-1; % precission
    r = randi(p,arr,1);
    r = (r-1) ./ (p-1); 
end 

