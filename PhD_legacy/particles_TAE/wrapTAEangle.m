 %wrap - Set radian angles in range [0,TAE_angle]
 

 function y = wrapTAEangle(x,TAE_angle)
 
 % Check number of parameters
 if nargin < 2,
   error('Incorrect number of parameters (type ''help <a href="matlab:help wrap">wrap</a>'' for details).');
 end
 

 % Determine angle in [0,TAE_angle]
 y = mod(x,TAE_angle);
 
