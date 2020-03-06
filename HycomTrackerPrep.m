function [V,G]=HycomTrackerPrep(varargin)
% HycomTrackerPrep Prep code for tracking particles through HYCOM model
% velocity fields. 
%
% Call as: [V,G]=HycomTrackerPrep('Url',url,'Level',level,'SubRegion',subregion,'ShiftLon',shiftlon)
%
% Inputs:
%   url       - (OPT) points to a 3D Hycom solution file on a THREDDS DATA SERVER:
%               (def='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015')
%               or to a file on local disk with u,v; use the NCSS THREDDS
%               service to downoad subsets of the full solutions.
%   level     - (OPT) is the vertical level to extract (def=1)
%   stride    - (OPT) spatial stride, basically the skip interval (def=1, no skip)
%   subregion - (OPT) is a vector with fields: 
%                    lon1  - lower-left longitude  
%                    lon2  - upper-right longitude 
%                    lat1  - lower-left latitude
%                    lat2  - upper-right latitude
%               (def=[-80  -50  5   35], Gulf oef Mexico)
%               OR: -1 to use the full grid in the netCDF file/URL
%   shiftlon  - (OPT) shift longitudes in HYCOM to -180:180. 
%               (def = true)
%      ndays  - number of days to access, from the start of the time in the
%               dataset at the url (def="all"
%   
% Outputs:
%   V - a struct V with fields, extracted from the url end-point:
%       u - east/west velocity field [nt X ne X nn]
%       v - north/south velocity field [nt X ne X nn]
%       t - time vector, in DATENUM format [nt]
% 
%   G -  a grid struct for fake finite element grid for drog2ddt
%       - x, y, z, e, bnd, etc... 
%
% Example:
%   [V,G]=HycomTrackerPrep;
%   IC=HycomTrackerIC;
%   R=HycomTracker(V,G,IC)
%

% Brian Blanton
% Renaissance Computing Institute
% The University of North Carolina at Chapel Hill
% Summer 2015
% Summer 2018

% Default propertyname values
url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
level=1;
%subregion=[lon1 lon2 lat1 lat2];
subregion =[-80  -50  5   35];
%subregion = -1;  % use all of the spatial region in the URL
PlotProgress=false;
%verbose=false;
stride=1;
shiftlon=true;
ndays=[];

if ~exist('ncgeodataset','file')
    error('This code requires nctoolbox for netCDF file access and processing. Install via GitHub from https://github.com/nctoolbox/nctoolbox')
end

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin)
  switch lower(varargin{k})
    case 'url'
      url=varargin{k+1};
      varargin([k k+1])=[];
    case 'level'
      level=varargin{k+1};
      varargin([k k+1])=[];
    case 'subregion'
      subregion=varargin{k+1};
      varargin([k k+1])=[];
    case 'plotprogress'
      subregion=varargin{k+1};
      varargin([k k+1])=[];
    case 'verbose'
      verbose=varargin{k+1};
      varargin([k k+1])=[];
    case 'stride'
      stride=varargin{k+1};
      varargin([k k+1])=[];
    case 'shiftlon'
      shiftlon=varargin{k+1};
      varargin([k k+1])=[];
    case 'ndays'
      ndays=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end     

% open the data pipe
fprintf('Contacting %s ...\n',url)
try 
    nc=ncgeodataset(url);
    %b=nc.variables{:}';
    %fprintf('Variables are: %s\n',b);
catch ME
     msg = sprintf('Failed to open netCDF link: %s.  Check link and server status.',url);
     causeException = MException('MATLAB:HycomTrackerPrep:ThreddsLink',msg);
     ME = addCause(ME,causeException);
     rethrow(ME)
end

fprintf('Got it.\nnc.location=%s\n',nc.location)

% get spatial variable objects
lon=nc{'Longitude'};
lat=nc{'Latitude'};
u=nc{'u'};
v=nc{'v'};
%ssh=nc{'ssh'}; 

lon_1d=cast(lon(1,:),'double');
lat_1d=cast(lat(:,1),'double');

% check lon var attributes for "modulo"
if strcmp(lon.attribute('modulo'),'360 degrees')
    shiftlon=true;
end

if shiftlon
    lon_1d=lon_1d-360;
end

% spatial subsetting region
if length(subregion)==4 
    lon1=subregion(1);
    lon2=subregion(2);
    lat1=subregion(3);
    lat2=subregion(4);
elseif subregion(1)==-1
    lon1=min(lon_1d);
    lon2=max(lon_1d);    
    lat1=min(lat_1d);
    lat2=max(lat_1d);
else
    error('Subregion option not valid.\n');
end

[~,ilon1]=min(abs(lon_1d-lon1));
[~,ilon2]=min(abs(lon_1d-lon2));
[~,ilat1]=min(abs(lat_1d-lat1));
[~,ilat2]=min(abs(lat_1d-lat2));

% create time vector
%datestr(nc.timeextent('u'))
Time=nc.gettimevar('MT');
units=Time.attribute('units');
base_date=datenum(units(12:end));
time=Time(:)+base_date;
dt=time(2)-time(1);
fprintf('%d time levels in %s\n',Time.size,url)

if isempty(ndays)
    ndays=length(time);
else
    if ndays > length(time)
        fprintf('Input ndays exceeds length of time at url. Setting to length(Time).\n');
        ndays=length(time);
    end
end

% if PlotProgress
%     figure 
%     fprintf('Extracting first time level of ssh ...\n');
%     sshd=ssh(1,:,:);
%     sshd=squeeze(cast(sshd,'double'));
%     pcolor(lon_1d,lat_1d,sshd)
%     shading flat
%     colorbar
%     axis('equal')
%     axis('tight')
%     caxis([-2 2])
%     colormap(parula(20))
%     line([lon1 lon2 lon2 lon1 lon1],[lat1 lat1 lat2 lat2 lat1])
%     title({url,['Full Domain, Sea Surface Height [m] @ ' datestr(time(1),2)]},'Interpreter','none')
%     print -dpng -r100 GlobalHycomGrid.png
%     
%     axis([lon1 lon2 lat1 lat2])
%     title({url,['Subregion Domain, Sea Surface Height [m] @ ' datestr(time(1),2)]},'Interpreter','none')
%     caxis([-1 1])
%     colormap(parula(10))
%     print -dpng -r100 SubRegionHycomGrid.png
% 
% end

% this is the subgrid
ilon=ilon1:stride:ilon2;
nlon=numel(ilon);
ilat=ilat1:stride:ilat2;
nlat=numel(ilat);
[lon_2d,lat_2d]=meshgrid(lon_1d(ilon),lat_1d(ilat));

% Create a fake ADCIRC grid for drog2ddt
fprintf('Creating fake finite element grid for HYCOM ...\n');
G.lon=lon_2d(:);
G.lat=lat_2d(:);
G.lo0=mean(G.lon);
G.la0=mean(G.lat);
G.lon1=lon1;
G.lon2=lon2;
G.lat1=lat1;
G.lat2=lat2;
G.z=NaN*ones(size(G.lat));
G.name='FakeHycomGrid';
G.e=elgen(nlon,nlat);
G.bnd=detbndy(G.e);
% G=belint(G);
% G=el_areas(G);
% project lon/lat to xy since u,v are in m/s
[G.x,G.y]=convll2m(G.lon,G.lat,G.lo0,G.la0);
G=belint(G);    % interpolation/basic functions
G=el_areas(G);  % element areas in m^2

% create velocity dataset for tracking
%fprintf('Extracting u,v at level=%d and ssh for %d time levels...\n',level,Time.size)
fprintf('Extracting u,v at level=%d for %d time levels...\n',level,Time.size)
ut=NaN*ones(Time.size,nlat,nlon);
vt=ut;
st=ut;
for it=1:ndays
    fprintf('Loading time level %d (%s)\n',it,datestr(time(it),2))
    ut(it,:,:)=squeeze(u(it,level,ilat,ilon)); 
    vt(it,:,:)=squeeze(v(it,level,ilat,ilon)); 
    %st(it,:,:)=squeeze(ssh(it,ilat,ilon)); 
end

%ut=cast(ut,'double');
%vt=cast(vt,'double');
%st=cast(st,'double');

V(Time.size).u=NaN;
V(Time.size).v=NaN;
V(Time.size).ssh=NaN;
V(Time.size).time=NaN;

fprintf('Restructuring arrays for tracker...\n')
for i=1:Time.size
    temp=squeeze(ut(i,:,:));
    V(i).u=temp(:);
    temp=squeeze(vt(i,:,:));
    V(i).v=temp(:);
    temp=squeeze(st(i,:,:));
    V(i).ssh=temp(:);
%    V(i).time=datetime(datevec(time(i)));
    V(i).time=time(i);
end

V(1).url=url;
V(1).PlotProgress=PlotProgress;


% plot first velocity time level
% if PlotProgress
%     hold on
%     quiver(lon_2d,lat_2d,squeeze(ut(1,:,:)),squeeze(vt(1,:,:)))
% end

