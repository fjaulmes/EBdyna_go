function s = makimaSlopes(delta)
% Derivative values for modified Akima cubic Hermite interpolation

% Akima's derivative estimate at grid node x(i) requires the four finite
% differences corresponding to the five grid nodes x(i-2:i+2).
%
% For boundary grid nodes x(1:2) and x(n-1:n), append finite differences
% which would correspond to x(-1:0) and x(n+1:n+2) by using the following
% uncentered difference formula correspondin to quadratic extrapolation
% using the quadratic polynomial defined by data at x(1:3)
% (section 2.3 in Akima's paper):
    n = numel(delta) + 1; % number of grid nodes x
    delta_0  = 2*delta(1)   - delta(2);
    delta_m1 = 2*delta_0    - delta(1);
    delta_n  = 2*delta(n-1) - delta(n-2);
    delta_n1 = 2*delta_n    - delta(n-1);
    delta = [delta_m1 delta_0 delta delta_n delta_n1];

% Akima's derivative estimate formula (equation (1) in the paper):
%
%       H. Akima, "A New Method of Interpolation and Smooth Curve Fitting
%       Based on Local Procedures", JACM, v. 17-4, p.589-602, 1970.
%
% s(i) = (|d(i+1)-d(i)| * d(i-1) + |d(i-1)-d(i-2)| * d(i))
%      / (|d(i+1)-d(i)|          + |d(i-1)-d(i-2)|)
%
% To eliminate overshoot and undershoot when the data is constant for more
% than two consecutive nodes, in MATLAB's 'makima' we modify Akima's
% formula by adding an additional averaging term in the weights.
% s(i) = ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) * d(i-1) +
%          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2) * d(i)  )
%      / ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) +
%          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2)
    weights = abs(diff(delta)) + abs((delta(1:end-1)+delta(2:end))/2);

    weights1 = weights(1:n);   % |d(i-1)-d(i-2)|
    weights2 = weights(3:end); % |d(i+1)-d(i)|
    delta1 = delta(2:n+1);     % d(i-1)
    delta2 = delta(3:n+2);     % d(i)

    weights12 = weights1 + weights2;
    s = (weights2./weights12) .* delta1 + (weights1./weights12) .* delta2;

% If the data is constant for more than four consecutive nodes, then the
% denominator is zero and the formula produces an unwanted NaN result.
% Replace this NaN with 0.
    s(weights12 == 0) = 0;
end