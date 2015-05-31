function IC=HycomTrackerIC(G)

% build out this function to suit needs.  The example below is a cluster
% of starting locations off the Texas coast:

lon1=264;
lon2=265;
lat1=25.5;
lat2=26.5;

dll=.1;
lon=lon1:dll:lon2;
lat=lat1:dll:lat2;
[IC.lon, IC.lat]=meshgrid(lon,lat);
IC.lon=IC.lon(:);
IC.lat=IC.lat(:);

% find initial elements
% fprintf('Locating initial positions in fake grid...\n')
% IC.j=findelem(G,IC.x,IC.y);