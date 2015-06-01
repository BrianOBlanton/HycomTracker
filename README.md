# HycomTracker
Demonstration MATLAB code for particle tracking with HYCOM model output posted to THREDDS Data Servers (TDS).  The default url is http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015. 


## nctoolbox
This code uses nctoolbox to access and extract variables from urls. You can get nctoolbox from https://github.com/nctoolbox/nctoolbox.  Either download as a zip file or clone in Desktop.  Then, add the path to MATLAB and run the setup code, something like: 

<pre>
addpath('<path_to_nctoolbox>/nctoolbox')
setup_nctoolbox; 
</pre>


There are 3 main codes: 

+ HycomTrackerPrep - Builds velocity arrays and fake HYCOM grid.  Default subregion is Caribbean and Gulf of Mexico, surface velocities. The default HYCOM URL is http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015

+ HycomTrackerIC - Generates initial particle locations.  Default is a square patch off the TX/MEX border.  User can edit/enhance as needed.

+ HycomTracker - This converts lon/lat coords to cartesian, calls drog2ddt, and inverts the projected particle locations back to lon/lat.

There are other supporting codes that handle grid structures and coordinate projections.  

## Run in MATLAB.


<pre>
% First step: generate the velocity arrays from HYCOM and fake HYCOM/finite element grid.  
% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015
[V,G]=HycomTrackerPrep;

% Then, build initial condition/locations:

IC=HycomTrackerIC;

% Pass to drog2ddt, through HycomTracker handler:
R=HycomTracker(V,G,IC);

% Then, make plots of R.lon and R.lat, etc...
% Subtract 360 from longitude to get it in the range -180->180; note the transpose.
url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
plot(R.lon'-360,R.lat')       % plot the trajectories
axis('equal')
line(R.lon(:,1)-360,R.lat(:,1),'Marker','.','Color','k','LineStyle','none')    % plot the initial positions
% plot some coastline data, if you have it.  If you have the mapping toolbox, then:
c=load('coast');   % note that this data is already in the range -180->180
line(c.long,c.lat,'Color','k')
axis([-99 -75 16 34])
title({'Trajectories in HYCOM Surface Velocity',url,[datestr(V(1).time,2) ' thru ' datestr(V(end).time,2) ]},'Interpreter','none')
</pre>

Here's an example plot:

![TestImage1](ExamplePlot.png "Example")
