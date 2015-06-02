function R=HycomTracker(V,G,IC)
% Demonstration MATLAB code for particle tracking with HYCOM model output 
% posted to THREDDS Data Servers (TDS).  The default url is:
%      http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015. 
%
% This code uses nctoolbox to access and extract variables from urls. 
% You can get nctoolbox from https://github.com/nctoolbox/nctoolbox.  
% Either download as a zip file or clone in Desktop.  Then, add the path 
% to MATLAB and run the setup code, something like: 
%     addpath('<path_to_nctoolbox>/nctoolbox')
%     setup_nctoolbox; 
%
% There are 3 main codes: 
% + HycomTrackerPrep - Builds velocity arrays and fake HYCOM grid.  
%   Default subregion is Caribbean and Gulf of Mexico, surface velocities. 
%   The default HYCOM URL is http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015
% + HycomTrackerIC - Generates initial particle locations.  
%   Default is a square patch off the TX/MEX border.  User can edit/enhance as needed.
% + HycomTracker - This converts lon/lat coords to cartesian, calls drog2ddt, 
%   and inverts the projected particle locations back to lon/lat.
%
%There are other supporting codes that handle grid structures and coordinate projections.
%
% Run in MATLAB.
%   First step: generate the velocity arrays from HYCOM and fake HYCOM/finite element grid.  
%      [V,G]=HycomTrackerPrep;
% Then, build initial condition/locations:
%      IC=HycomTrackerIC;
% Then, pass to drog2ddt, through HycomTracker handler:
%      R=HycomTracker(V,G,IC);
%


options.draw=true;
t1d=V(1).time;
t2d=V(end).time;
idt=1;
dt=24;  % time step in hours


% project lon/lat to xy since u,v are in m/s
[G.x,G.y]=convll2m(G.lon,G.lat,G.lo0,G.la0);
G=belint(G);    % interpolation/basic functions
G=el_areas(G);  % element areas in m^2
[IC.x,IC.y]=convll2m(IC.lon,IC.lat,G.lo0,G.la0);




figure 
line(IC.x-360,IC.y,'Marker','.','Color','k','LineStyle','none')    
axis([min(G.x) max(G.x) min(G.y) max(G.y)])

[R.xx,R.yy,R.tt,R.uu,R.vv]=drog2ddt(G,t1d,t2d,dt,idt,IC.x,IC.y,V,options);

line(R.xx',R.yy')

% invert the forward projection for drogue locations
[R.lon,R.lat]=convm2ll(R.xx,R.yy,G.lo0,G.la0);

%%% Then, make plots of R.lon and R.lat, etc...
figure
% plot the trajectories
% Subtract 360 from longitude to get it in the range -180->180; note the transpose.
plot(R.lon'-360,R.lat')      
axis('equal')

% plot the initial positions
line(R.lon(:,1)-360,R.lat(:,1),'Marker','.','Color','k','LineStyle','none')    
line(R.lon(:,end)-360,R.lat(:,end),'Marker','.','Color','r','LineStyle','none','MarkerSize',14)    

% plot some coastline data, if you have it.  If you have the mapping toolbox, then:
if exist('coast.mat','file')
   c=load('coast');      
   % note that this data is already in the range -180->180
   line(c.long,c.lat,'Color','k')
end

axis([-99 -75 16 34])
title({'Trajectories in HYCOM Surface Velocity',...
        V(1).url,...
        [datestr(V(1).time,2) ' thru ' datestr(V(end).time,2) ]},'Interpreter','none')
