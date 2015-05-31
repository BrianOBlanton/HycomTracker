function [x,y]=convll2m(lo,la,lo0,la0) 
%CONVLL2M convert long/lat (deg) to cartesian coordinates (meters)
% CONVLL2M converts from long/lat coordinates (in degrees) to
% meters.  
%
%    Input: lo  - longitudes in degrees
%           la  - latitudes in degrees
%           lo0 - reference longitude
%           la0 - reference latitude
%  Outputs: x   - east coordinate in meters
%           y   - west coordinate in meters
%
%  Call as:  [x,y]=convll2m(lo,la,lo0,la0);
%
% Written by : Brian O. Blanton
%              Summer 2000

 
if nargin < 4
   error('Too few arguments to CONVLL2M')
end
if nargout ~=2
   error('CONVLL2M requires two return arguments.')
end

R=6367500;  % Approx earth radius in meters

% Check sizes of lo,la
if ~all(size(lo)==size(la))
   error('Sizes of lo and la are NOT equal')
end

deg2rad=pi/180;
fac=R*cos(la0*deg2rad);


x=fac*lo*deg2rad;
y=fac*log((1+sin(la*deg2rad))./(cos(la*deg2rad)));

return

%        Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
%        Summer 2000
%


