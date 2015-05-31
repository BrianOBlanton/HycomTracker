function R=HycomTracker(V,G,IC)

options.draw=true;
t1d=V(1).time;
t2d=V(end).time;
idt=1;
dt=24;  % time step in hours


% project lon/lat to xy since u,v are in m/s
G.lo0=mean(G.lon);
G.la0=mean(G.lat);
[G.x,G.y]=convll2m(G.lon,G.lat,G.lo0,G.la0);
G=belint(G);
G=el_areas(G);

[IC.x,IC.y]=convll2m(IC.lon,IC.lat,G.lo0,G.la0);

[R.xx,R.yy,R.tt,R.uu,R.vv]=drog2ddt(G,t1d,t2d,dt,idt,IC.x,IC.y,V,options);

% invert the forward projetion for drogue locations
[R.lon,R.lat]=convm2ll(R.xx,R.yy,G.lo0,G.la0);

