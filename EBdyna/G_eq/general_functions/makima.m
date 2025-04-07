function vq = makima(x,v,xq)
%MAKIMA   Modified Akima piecewise cubic Hermite interpolation.
%
%       We modify Akima's formula to eliminate overshoot and undershoot
%       when the data is constant for more than two consecutive nodes.
%
%   vq = MAKIMA(x,v,xq) interpolates the data (x,v) at the query points xq
%   using a modified Akima cubic Hermite interpolation formula.
%
%   References:
%
%       H. Akima, "A New Method of Interpolation and Smooth Curve Fitting
%       Based on Local Procedures", JACM, v. 17-4, p.589-602, 1970.
%
%   MAKIMA vs. PCHIP vs. SPLINE:
%
%     - MAKIMA is a middle ground between SPLINE and PCHIP:
%       It has lower-amplitude wiggles than SPLINE, but is not as agressive
%       at reducing the wiggles as PCHIP.
%     - MAKIMA and SPLINE generalize to n-D grids.
%
%   Example: No overshoot and undershoot with MAKIMA when the data is
%            constant for more than two consecutive nodes
%
%       x = 1:7;
%       v = [-1 -1 -1 0 1 1 1];
%       xq = 0.75:0.05:7.25;
%       vqs = spline(x,v,xq);
%       vqp = pchip(x,v,xq);
%       vqa = akima(x,v,xq);
%       vqm = makima(x,v,xq);
%
%       figure
%       plot(x,v,'ko','LineWidth',2,'MarkerSize',10), hold on
%       plot(xq,vqp,'LineWidth',4)
%       plot(xq,vqs,xq,vqa,xq,vqm,'LineWidth',2)
%       legend('Data','''pchip''','''spline''','Akima''s formula',...
%           '''makima''','Location','SE')
%       title('''makima'' has no overshoot in x(1:3) and x(5:7)')
%
%   Example: Compare MAKIMA with AKIMA, SPLINE, and PCHIP
%
%       x = [1 2 3 4 5 5.5 7 8 9 9.5 10];
%       v = [0 0 0 0.5 0.4 1.2 1.2 0.1 0 0.3 0.6];
%       xq = 0.75:0.05:10.25;
%       vqs = spline(x,v,xq);
%       vqp = pchip(x,v,xq);
%       vqa = akima(x,v,xq);
%       vqm = makima(x,v,xq);
%
%       figure
%       plot(x,v,'ko','LineWidth',2,'MarkerSize',10), hold on
%       plot(xq,vqp,'LineWidth',4)
%       plot(xq,vqs,xq,vqa,xq,vqm,'--','LineWidth',2)
%       legend('Data','''pchip''','''spline''','Akima''s formula','''makima''')
%       title('Cubic Hermite interpolation in MATLAB')
%
%   See also AKIMA, SPLINE, PCHIP.


%   mailto: cosmin.ionita@mathworks.com

    assert(isvector(x) && isvector(v) && (numel(x) == numel(v)))
    assert(numel(x) > 2)
    x = x(:).'; v = v(:).'; % Same shapes as in pchip.m and spline.m

% Compute modified Akima slopes
    h = diff(x);
    delta = diff(v)./h;
    slopes = makimaSlopes(delta);

% Evaluate the piecewise cubic Hermite interpolant
    vq = ppval(pwch(x,v,slopes,h,delta),xq);

end