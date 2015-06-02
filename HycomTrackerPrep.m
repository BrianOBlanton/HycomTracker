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
subregion =[262  295     7   32];
PlotProgress=true;
verbose=false;

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
      subregion=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;      
     
% open the data pipe
try 
    url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1/2015';
    nc=ncgeodataset(url);
    fprintf('nc.location=%s\n',nc.location)
    fprintf('nc.netcdf=%s\n',nc.netcdf)
    %b=nc.variables{:}';
catch
    fprintf('Failed to open netCDF link: %s.  Check link and server status.\n',url)
end

% get spatial variables
lon=nc{'Longitude'};
lat=nc{'Latitude'};
u=nc{'u'};
v=nc{'v'};
ssh=nc{'ssh'}; 
 
lon_1d=cast(lon(1,:),'double');
lat_1d=cast(lat(:,1),'double');

% define spatial subsetting region
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
lon_2d=cast(lon(ilat1:ilat2,ilon1:ilon2),'double');
lat_2d=cast(lat(ilat1:ilat2,ilon1:ilon2),'double');
[nvert,nhoriz]=size(lon_2d);

% Create a fake ADCIRC grid for drog2ddt
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
it=1:Time.size;
%it=1:2;
ut=squeeze(u(it,level,ilat1:ilat2,ilon1:ilon2)); 
vt=squeeze(v(it,level,ilat1:ilat2,ilon1:ilon2)); 
st=squeeze(ssh(it,level,ilat1:ilat2,ilon1:ilon2)); 
ut=cast(ut,'double');
vt=cast(vt,'double');
st=cast(st,'double');

for i=1:length(it)
    temp=squeeze(ut(i,:,:));
    V(i).u=temp(:);
    temp=squeeze(vt(i,:,:));
    V(i).v=temp(:);
    temp=squeeze(st(i,:,:));
    V(i).st=temp(:);
    V(i).time=time(i);
end

V(1).url=url;
V(1).PlotProgress=PlotProgress;


% plot first velocity time level
if PlotProgress
    hold on
    quiver(lon_2d,lat_2d,ut,vt)
end

