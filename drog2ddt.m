function [xx,yy,tt,uu,vv]=drog2ddt(fem_grid_struct,t1d,t2d,dt,idt,xi,yi,V,options)
%DROG2DDT track drogues in a 2-D FEM domain, time-stepping version
% DROG2DDT tracks particles through a discrete squence of
% 2-D velocity fields, vertically averaged for example.  The
% integrator is a 2nd order Runge-Kutta (mid-point) method.
%
% Inputs: fem_grid_struct - FEM domain structure (from LOADGRID)
%                           Horizontal coordinates are in meters (CARTESIAN).
%         t1,t2           - integration end-points; both t1 & t2
%                           must lie within the min and max time
%                           in the velocity sequence.  The time units
%                           are days!!  (as in MATLAB's datenum
%                           function, for example).
%         dt              - integration time step; this need NOT
%                           match the velocity sequence step.  The units
%                           are specified in HOURS!!
%         idt             - output interval; this is the
%                           "frequency" with which to store
%                           computed step.  idt=1 means store 
%                           each step; idt=2 means store every
%                           other step, etc.
%         xi,yi           - initial drogue positions
%         V               - the sequence of velocity fields. 
%                           Type "drog2ddt('velhelp')" for information.
%                           Units must be in meters/sec.  The time field
%                           of the velocity sequence must be in days!!,
%                           consistent with the t1,t2 tracking interval.
%         options         - "options" structure; valid fields are:
%                           "draw" - place updated locations on figure
%
% Outputs: xx,yy - arrays of drogue trajectories.  The size of the
%                  arrays is the size of xi,yi, plus an extra 
%                  dimension for the time history.  This dimension
%                  is the number of iterations used, possibly 
%                  subsampled at idt intervals.  Example: if
%                  xi is 10x10, and t1=0,t2=10,dt=1,idt=2, then the 
%                  size of xx is 10x10x6.
%
%          tt    - times at which outputs are stored
%
%          uu,vv - along-track velocity in the same format as xx,yy
% 
% Call as: [xx,yy,tt,uu,vv]=drog2ddt(fem_grid_struct,t1,t2,dt,idt,xi,yi,V,options);
% 

% Written by: Brian Blanton, Fall 99
%             added output of along-track u,v 31 Jan 03
if nargin==0
   disp('Call as:  [xx,yy]=drog2ddt(fem_grid_struct,t1,t2,dt,idt,xi,yi,V,options)')
   return
end

if nargin==1&strcmp(fem_grid_struct,'velhelp')
   velhelp('Help');
   return
end

if nargin<8
   error('Too few arguments.')
end

% First argument to drog2ddt must be a grid structure
if ~isstruct(fem_grid_struct)
   error('First argument to DROG2DDT must be a fem_grid_struct')
end
if ~is_valid_struct(fem_grid_struct)
   error('fem_grid_struct to DROG2DDT NOT valid.')
end


% Check options structure
if ~exist('options')
   options.draw=0;
else
   if ~isstruct(options)
     error('Options argument to DROG2DDT must be a structure')
   end
end

% Verify velocity array of structures
if ~isstruct(V)
   velhelp('Error');
   error('Velocity array must be a "struct array".')
end
%Check that field names of V include atleast u,v,time
fnames=fieldnames(V);
if ~all(any([strcmp(fnames,'u') strcmp(fnames,'v') strcmp(fnames,'time')]))
   error(['Velocity array lacks a needed field name.'])
end
%Loop over number of velocity snapshots and verify
sgx=size(fem_grid_struct.x);
for i=1:length(V)
   if ~all(size(V(i).u)==sgx)
      error(['u field for V(' int2str(i) ') not same size as grid.'])
   end
   if ~all(size(V(i).v)==sgx)
      error(['v field for V(' int2str(i) ') not same size as grid.'])
   end
end


% convert times from days to hours
t1=t1d*24;
t2=t2d*24;
for i=1:length(V)
   V(i).time=V(i).time*24;
end


% check timestep
if dt<eps
   error('Timestep<eps???  You''re kidding, right???')
end

disp(['Tracking Start: '  datestr(t1d,0) ])
disp(['Tracking End  : '  datestr(t2d,0) ])

% Error-check the input times relative to the timestamps
% in the velocity cell array V.  It is ASSUMED that the 
% velocity slices are in temporal order within the cells
timemin=V(1).time;
timemax=V(length(V)).time;
if t1 < timemin
   error('T1 less than minimum time in velocity arrays.')
elseif t2 > timemax
   error('T2 exceeds maximum time in velocity arrays.')
end


% Check sizes of input drogue positions
if ~all(size(xi)==size(yi))
   error('Sizes of initial drog position arrays must be equal.')
end
[mdrog,ndrog]=size(xi);
xi=xi(:);yi=yi(:);

% Extract a time vector for the times at which the velocity
% slices are available.
timevec=[V.time];

% Attach element finding arrays to fem_grid_struct
if ~isfield(fem_grid_struct,'A') |~isfield(fem_grid_struct,'B')|...
   ~isfield(fem_grid_struct,'A0')|~isfield(fem_grid_struct,'T')
   fem_grid_struct=belint(fem_grid_struct);
   disp('   BELINT info added to fem_grid_struct')
end
if ~isfield(fem_grid_struct,'ar')
   disp('   EL_AREAS info added to fem_grid_struct')
   fem_grid_struct=el_areas(fem_grid_struct);
end

% Locate initial positions in grid
% j will be the array that keeps track of the
% current element number of each drog.  A NaN will
% be used to indicate drog in-activity, either because
% the drog is initially out-of-bounds, or because the
% drogue has left the domain during tracking.
j=findelem(fem_grid_struct,[xi yi]);

% Allocate space for time history of positions
tt=t1:dt:t2;tt=tt(1:idt:length(tt));
xx=NaN*(ones(size(tt))'*ones(size(xi')))';
yy=xx;
uu=xx;
vv=xx;

% 
xx(:,1)=xi;
yy(:,1)=yi;
[ut,vt]=vel_interp(fem_grid_struct,xi,yi,j,V,timevec,timevec(1));
uu(:,1)=ut;
vv(:,1)=vt;

%Draw initial positions on screen
if options.draw
   hdrog=line(xi,yi,'LineStyle','none','Marker','.', ...
                    'MarkerSize',14,'Color','r');
   drawnow
end

% Loop over time;
time=t1;
iter=1;
ic=0;
dtsecs=dt*3600;
disp(['Starting: [' int2str(iter) ' ' datestr(time/24,0) ']'])
xnew=xi;
ynew=yi;


while time<t2
   disp(['   Integrating: [' int2str(iter) ' ' datestr(time/24,0) ']'])
   % 
%   xnow=xx(:,iter);
%   ynow=yy(:,iter);
   xnow=xnew;
   ynow=ynew;
   
%   xnew=xnow;   % Propagate previously known positions, in case a 
%   ynew=ynow;   % drogue is eliminated.

   igood=find(~isnan(j));
   % If j contains NaN's, these drogues are have exited the
   % domain at some previous timestep. 
   
   if isempty(igood)
      disp('All drogues eliminated')
      disp(['Ending: [' int2str(iter) ' ' datestr(time/24,0) ']'])
      return
   end
   
   % Extract drogues currently in domain
   jgood=j(igood);xgood=xnow(igood);ygood=ynow(igood);
   
   [xnext,ynext,jnext]=track(fem_grid_struct,jgood,xgood,ygood,V,timevec,time,dtsecs);
   j(igood)=jnext;
   xnew(igood)=xnext;
   ynew(igood)=ynext;

   time=time+dt;
   iter=iter+1;   

   % Store this timestep if its time
   if rem(iter,idt)==0
   %keyboard
      ic=ic+1;
      [ut,vt]=vel_interp(fem_grid_struct,xnext,ynext,jnext,V,timevec,time);
      xx(:,ic+1) =xnew;
      yy(:,ic+1) =ynew;
      uu(igood,ic+1)  =ut;
      vv(igood,ic+1)  =vt;
      tt(ic+1)    =time;
   end
   
   %Update positions on screen
   if options.draw
      delete(hdrog)
      hdrog=line(xnew,ynew,'LineStyle','none','Marker','.', ...
	   'MarkerSize',14,'Color','r');
      drawnow
   end
   
end

% Prepare for return
disp(['Ending: [' int2str(iter) ' ' datestr(time/24,0) ']'])
if mdrog==1
   xx=reshape(xx,ndrog,ic+1);
   yy=reshape(yy,ndrog,ic+1);
   uu=reshape(uu,ndrog,ic+1);
   vv=reshape(vv,ndrog,ic+1);
elseif ndrog==1
   xx=reshape(xx,mdrog,ic+1);
   yy=reshape(yy,mdrog,ic+1);
   uu=reshape(uu,mdrog,ic+1);
   vv=reshape(vv,mdrog,ic+1);
else
   xx=reshape(xx,mdrog,ndrog,ic+1);
   yy=reshape(yy,mdrog,ndrog,ic+1);
   uu=reshape(uu,mdrog,ndrog,ic+1);
   vv=reshape(vv,mdrog,ndrog,ic+1);
end

return



% PRIVATE FUNCTIONS

function jnew=locate_drog(fem_grid_struct,x,y,j)
jnew=j;
% See if drogs are still in Previously known element
idx=belel(fem_grid_struct,j,[x y]);
inotfound=find(idx==0);
% Get new elements, if not in previously known element
if ~isempty(inotfound)
   idx=inotfound;
   jnew(idx)=findelem(fem_grid_struct,[x(idx) y(idx)]);
end

function [xnew,ynew,jnew]=track(fem_grid_struct,j,x,y,V,timevec,t,dt)
jnew=j;
% Mid-Point (RK2) step
[u,v]=vel_interp(fem_grid_struct,x,y,j,V,timevec,t);
uk1=dt*u;vk1=dt*v;
xinter=x+.5*uk1;
yinter=y+.5*vk1;
jinter=locate_drog(fem_grid_struct,xinter,yinter,j);
[uinter,vinter]=vel_interp(fem_grid_struct,xinter,yinter,jinter,V,timevec,t+.5*dt/3600);
uk2=dt*uinter;
vk2=dt*vinter;
xnew=x+uk2;
ynew=y+vk2;
         
% If NaN is in j, then drog has left domain.  
% insert last known location into arrays
jnew=locate_drog(fem_grid_struct,xnew,ynew,j);
inan=find(isnan(jnew));;
if ~isempty(inan)
   xnew(inan)=x(inan);
   ynew(inan)=y(inan);
end



function retval=belel(fem_grid_struct,j,xylist)
%BELEL - determine if points are in elements
% BELEL 
tol=eps*10000000;
phi=basis2d(fem_grid_struct,xylist,j);
test=phi>=-tol & phi<=1+tol;
retval=all(test'==1);

function [u,v]=vel_interp(fem_grid_struct,x,y,j,V,timevec,t)
% Get the velocities at this time
% Temporal interpolation of velocity slices to this time.
it1=find(t<=timevec);it1=it1(1);
it2=find(t>=timevec);it2=it2(length(it2));
if it1==it2
   tfac1=1;tfac2=0;
else
   tfac1=(timevec(it1)-t)/(timevec(it1)-timevec(it2));
   tfac2=1-tfac1;
end
u1=V(it1).u;
u2=V(it2).u;
v1=V(it1).v;
v2=V(it2).v;

% Depending on the number of particles to track, the
% interpolation to (x,y,t) is done one of two ways.
% It is not ovvious, but the flop savings can be huge.
if length(j)>150
   % Interpolate in time first,...
   uu=tfac1*u2 + tfac2*u1;
   vv=tfac1*v2 + tfac2*v1;
   % ... then, space
   u=interp_scalar(fem_grid_struct,uu,x,y,j);
   v=interp_scalar(fem_grid_struct,vv,x,y,j);
else
   % Interpolate in space first, at the 2 time levels
   uu1=interp_scalar(fem_grid_struct,u1,x,y,j);
   vv1=interp_scalar(fem_grid_struct,v1,x,y,j);
   uu2=interp_scalar(fem_grid_struct,u2,x,y,j);
   vv2=interp_scalar(fem_grid_struct,v2,x,y,j);
   % Then, interpolate BETWEEN time levels
   u=tfac1*uu2 + tfac2*uu1;
   v=tfac1*vv2 + tfac2*vv1;
end


function velhelp(ftitle)
str={'This describes the format of the array which contains the velocity',...
'sequence.  The array V must be an array of structures, one',...
'structure per timestep, with the following format;',...
'   V(1).u',...
'   V(1).v',...
'   V(1).time',...
'   ...',...
'   ...',...',...
'   V(n).u',...
'   V(n).v',...
'   V(n).time',...
'where n is the number of discrete velocity fields in V.  ',...
'The fields u,v,time must be part of each structure, and ',...
'the names of the fields MUST be u,v,time.  ',...
' ',...
'The fields u,v must be the same size as the x,y,z',...
'coordinate arrays in the fem_grid_struct.  The time',...
'field is the timestamp for the velocity snapshot in HOURS.'};
h=helpdlg(str,ftitle);



%
%        Brian O. Blanton
%        Department of Marine Sciences
%        Ocean Processes Numerical Modeling Laboratory
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
%        September  1999
%
