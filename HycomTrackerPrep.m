function [V,G]=HycomTrackerPrep(varargin)
%
% Call as: [V,G]=HycomTrackerPrep(url,'Level',ilevel,'SubRegion',subregion)
%
% where:
%   url (OPT) points to a 3D Hycom solution file on a THREDDS DATA SERVER:
%      Default ='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
%   ilevel (OPT) is the vertical level to extract (def=1)
%   subregion (OPT) is a vector with fields: 
%       lon1 - lower-left longitude
%       lon2 - upper-right longitude
%       lat1 - lower-left latitude
%       lat2 - upper-right latitude
%
% Returns a struct V with fields, extracted from the url end-point:
%   u - east/west velocity field [nt X ne X nn]
%   v - north/south velocity field [nt X ne X nn]
%   t - time vector, in DATENUM format [nt]
% 
% Also returns G, a grid struct for fake finite element grid for drog2ddt
%       - x, y, z, e, bnd, etc... 
%
% Pass both V and G to HycomTrackerIC to generate initial conditions for
% particle locations.
%

% Default propertyname values
url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
level=1;
%subregion=[lon1 lon2 lat1 lat2];
subregion=[262 295 7 32];

%verbose=true;

if ~exist('ncgeodataset')
    error('This code requires nctoolbox for netCDF file access and processing. Install via GitHub from https://github.com/nctoolbox/nctoolbox')
end

try 
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
    nc=ncgeodataset(url);
catch
    fprintf('Failed to open netCDF link: %s.  Check link and server status.\n',url)
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
    otherwise
      k=k+2;
  end;
end;      
     
% get spatial variables
lon=nc{'Longitude'};
lat=nc{'Latitude'};
u=nc{'u'};
v=nc{'v'};

lon_1d=cast(lon(1,:),'double');
lat_1d=cast(lat(:,1),'double');

% define spatial subsetting region
lon1=subregion(1);
lon2=subregion(2);
lat1=subregion(3);
lat2=subregion(4);

[~,ilon1]=min(abs(lon_1d-lon1));
[~,ilon2]=min(abs(lon_1d-lon2));
[~,ilat1]=min(abs(lat_1d-lat1));
[~,ilat2]=min(abs(lat_1d-lat2));

% this is the subgrid
lon_2d=cast(lon(ilat1:ilat2,ilon1:ilon2),'double');
lat_2d=cast(lat(ilat1:ilat2,ilon1:ilon2),'double');
[nvert,nhoriz]=size(lon_2d);

% Create a fake ADCIRC grid for drog2ddt
G.lon=lon_2d(:);
G.lat=lat_2d(:);
G.z=NaN*ones(size(G.lat));
G.name='FakeHycomGrid';
G.e=elgen(nhoriz,nvert);
G.bnd=detbndy(G.e);
% G=belint(G);
% G=el_areas(G);


% create velocity dataset for tracking
ut=squeeze(u(1,level,ilat1:ilat2,ilon1:ilon2)); 
vt=squeeze(v(1,level,ilat1:ilat2,ilon1:ilon2)); 
ut=cast(ut,'double');
vt=cast(vt,'double');

% create time vector
%datestr(nc.timeextent('u'))
Time=nc.gettimevar('MT');
units=Time.attribute('units');
base_date=datenum(units(12:end));
time=Time(:)+base_date;

for i=1:length(time)
    temp=squeeze(ut(i,:,:));
    V(i).u=temp(:);
    temp=squeeze(vt(i,:,:));
    V(i).v=temp(:);
    V(i).time=time(i);
end



