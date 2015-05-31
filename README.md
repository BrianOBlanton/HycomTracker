# HycomTracker
Demonstration MATLAB code for particle tracking with HYCOM model output

There are 3 main codes: 

1) HycomTrackerPrep - Builds velocity arrays and fake HYCOM grid.  Default subregion is Caribbean and Gulf of Mexico, surface velocities. The default HYCOM URL is http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015

2) HycomTrackerIC - Generates initial particle locations.  Default is a square patch off the TX/MEX border.  User can edit/enhance as needed.

3) HycomTracker - This converts lon/lat coords to cartesian, calls drog2ddt, and inverts the projected particle locations back to lon/lat.


## Run in MATLAB.
First step: generate the velocity arrays from HYCOM and fake HYCOM/finite element grid.  
http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015

>> [V,G]=HycomTrackerPrep;


>> IC=HycomTrackerIC;

>> R=HycomTracker(V,G,IC);

Then, make plots of R.lon and R.lat, etc...


