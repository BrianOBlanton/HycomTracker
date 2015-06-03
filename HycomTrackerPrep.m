function [V,G]=HycomTrackerPrep(varargin)
% HycomTrackerPrep Prep code for tracking particles through HYCOM model
% velocity fields.
%
% Call as: [V,G]=HycomTrackerPrep('Url',url,'Level',level,'SubRegion',subregion)
%
% Inputs:
%   url       - (OPT) points to a 3D Hycom solution file on a THREDDS DATA SERVER:
%               (def='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015')
%   level     - (OPT) is the vertical level to extract (def=1)
%   stride    - (OPT) spatial stride, basically the skip interval (def=1, no skip)
%   subregion - (OPT) is a vector with fields: 
%                    lon1  - lower-left longitude
%                    lon2  - upper-right longitude
%                    lat1  - lower-left latitude
%                    lat2  - upper-right latitude
%               (def=[262  295  7   32], Gulf oef Mexico)
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

% Default propertyname values
url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
level=1;
%subregion=[lon1 lon2 lat1 lat2];
subregion =[262  295     7   32];
PlotProgress=false;
verbose=false;
stride=1;

if ~exist('ncgeodataset','file')
    error('This code requires nctoolbox for netCDF file access and processing. Install via GitHub from https://github.com/nctoolbox/nctoolbox')
end

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'url',
      url=varargin{k+1};
      varargin([k k+1])=[];
    case 'level',
      level=varargin{k+1};
      varargin([k k+1])=[];
    case 'subregion',
      subregion=varargin{k+1};
      varargin([k k+1])=[];
    case 'plotprogress',
      subregion=varargin{k+1};
      varargin([k k+1])=[];
    case 'verbose',
      verbose=varargin{k+1};
      varargin([k k+1])=[];
    case 'stride',
      stride=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;      
     
% open the data pipe
fprintf('Trying to contact %s ...\n',url)
try 
    nc=ncgeodataset(url);
    %b=nc.variables{:}';
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
ssh=nc{'ssh'}; 

lon_1d=cast(lon(1,:),'double');
lat_1d=cast(lat(:,1),'double');

% spatial subsetting region
lon1=subregion(1);
lon2=subregion(2);
lat1=subregion(3);
lat2=subregion(4);

% create time vector
%datestr(nc.timeextent('u'))
Time=nc.gettimevar('MT');
units=Time.attribute('units');
base_date=datenum(units(12:end));
time=Time(:)+base_date;
fprintf('%d time levels in %s\n',Time.size,url)

if PlotProgress
    figure 
    fprintf('Extracting first time level of ssh ...\n');
    sshd=ssh(1,:,:);
    sshd=squeeze(cast(sshd,'double'));
    pcolor(lon_1d,lat_1d,sshd)
    shading flat
    colorbar
    axis('equal')
    axis('tight')
    caxis([-2 2])
    colormap(parula(20))
    line([lon1 lon2 lon2 lon1 lon1],[lat1 lat1 lat2 lat2 lat1])
    title({url,['Full Domain, Sea Surface Height [m] @ ' datestr(time(1),2)]},'Interpreter','none')
    print -dpng -r100 GlobalHycomGrid.png
    
    axis([lon1 lon2 lat1 lat2])
    title({url,['Subregion Domain, Sea Surface Height [m] @ ' datestr(time(1),2)]},'Interpreter','none')
    caxis([-1 1])
    colormap(parula(10))
    print -dpng -r100 SubRegionHycomGrid.png

end

[~,ilon1]=min(abs(lon_1d-lon1));
[~,ilon2]=min(abs(lon_1d-lon2));
[~,ilat1]=min(abs(lat_1d-lat1));
[~,ilat2]=min(abs(lat_1d-lat2));

% this is the subgrid
ilon=ilon1:stride:ilon2;
ilat=ilat1:stride:ilat2;
lon_2d=cast(lon(ilat,ilon),'double');
lat_2d=cast(lat(ilat,ilon),'double');
[nvert,nhoriz]=size(lon_2d);

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
G.e=elgen(nhoriz,nvert);
G.bnd=detbndy(G.e);
% G=belint(G);
% G=el_areas(G);

% create velocity dataset for tracking
fprintf('Extracting u,v at level=%d and ssh for %d time levels...\n',level,Time.size)
ut=NaN*ones(Time.size,nvert,nhoriz);
vt=ut;
st=ut;
for it=1:Time.size
    fprintf('Loading time level %d (%s)\n',it,datestr(time(it),2))
    ut(it,:,:)=squeeze(u(it,level,ilat,ilon)); 
    vt(it,:,:)=squeeze(v(it,level,ilat,ilon)); 
    st(it,:,:)=squeeze(ssh(it,ilat,ilon)); 
end

ut=cast(ut,'double');
vt=cast(vt,'double');
st=cast(st,'double');

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
    V(i).time=time(i);
end

V(1).url=url;
V(1).PlotProgress=PlotProgress;

% plot first velocity time level
if PlotProgress
    hold on
    quiver(lon_2d,lat_2d,squeeze(ut(1,:,:)),squeeze(vt(1,:,:)))
end

