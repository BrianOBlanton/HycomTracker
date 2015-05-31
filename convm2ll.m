function [lo,la]=convm2ll(x,y,lo0,la0)
%CONVM2LL convert cartesian coordinates (meters) to long/lat (deg)
% CONVM2LL converts from cartesian coordinates in meters to
% long/lat coordinates by inverting the conversion formulas 
% used in the routine CONVLL2M.  This is the inverse of
% converting long/lat coordinates to cartesian coordinates in
% meters. The x,y  coordinates MUST have been computed using
% CONVLL2M, and the  reference long/lat pair MUST be the same
% as those used in the  "forward" conversion".  Any
% translation in the x,y coordinates that might have been
% done MUST be "undone" prior to conversion back to long/lat
% coordinates.  The inversion for latitude is done with 5 
% Newton-Raphson iterations with an initial guess of the
% reference latitude.  The longitude is trivial.
%
%  Inputs: x,y - cartesian coordinates in meters
%          lo0,la0 - reference long/lat values
%
% Outputs: lo,la - long/lat coordinates.
%
% Call as:  [lo,la]=convm2ll(x,y,lo0,la0);

if nargin < 4
   error('Too few arguments to CONVM2LL')
end
if nargout ~=2
   error('CONVM2LL requires two return arguments.')
end

% Check sizes of lo,la
if ~all(size(x)==size(y))
   error('Sizes of lo and la are NOT equal')
end

R=6367500;  % Approx earth radius in meters
deg2rad=pi/180;
fac=R*cos(la0*deg2rad);

% Invert for lo (easy)
lo=x/fac/deg2rad;

% Invert for la (Newton iteration for nonlinear eqn)
% The eqn is y=fac*log((1+sin(la))./(cos(la)));
% where la is in radians and fac is 
%        R*cos(la0)
% where R is the earth's radius in meters and la0
% is the reference latitude.
% Rearange this into 
%                           1+sin(la)
%   f(la) = exp(y/fac) - --------------
%                            cos(la)
% and find la where f(la)==0
%
%                1+sin(la)
% f'(la) = -1 - ----------- * tan(la)
%                 cos(la)
%
phin=la0*ones(size(x))*deg2rad;
CC=exp(y/fac);
% 5 Newton-Raphson iterations, assume convergence!!
for i=1:5
   phinm1=phin;
   numer=CC - (1+sin(phinm1))./(cos(phinm1));
   denom= -1- (1+sin(phinm1)).*tan(phinm1)./cos(phinm1);
   phin=phinm1-numer./denom;
end

la=phin/deg2rad;



