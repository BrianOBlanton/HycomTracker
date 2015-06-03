function IC=HycomTrackerIC
%HycomTrackerIC - Generates initial particle locations for HycomTracker
% This is a completely user-defined code.  The returned struct needs to
% have two fields:, lon and lat as equally sized column vectors. It does
% not matter how the locations are generated as long as the coordinates are
% returned as 1-d vectors.  The drog2ddt code will determine if the
% locations are within the subregion grid.  Particles that are "on land"
% but within the subregion grid will not be tracked.
%
% Example:
%   [V,G]=HycomTrackerPrep;
%   IC=HycomTrackerIC;
%   R=HycomTracker(V,G,IC)

% build out this function to suit needs.  The example below is a 
% 1deg X 1deg cluster of starting locations off the Texas coast:

lon1=264;
lon2=265;
lat1=25.5;
lat2=26.5;

% separation in lon/lat
dll=.1;

lon=lon1:dll:lon2;
lat=lat1:dll:lat2;
[IC.lon, IC.lat]=meshgrid(lon,lat);
IC.lon=IC.lon(:);
IC.lat=IC.lat(:);

% find initial elements
% fprintf('Locating initial positions in fake grid...\n')
% IC.j=findelem(G,IC.x,IC.y);