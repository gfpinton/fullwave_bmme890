function [dY dZ] = launchTotalFullWaveRbuild2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap,rhomap,Amap,boveramap,incoords,outcoords,icmat)

nbdy = 40; % number of boundary points for PML

lambda = c0/omega0*2*pi;
lambda = c0/omega0*2*pi;
nY = round(wY/lambda*ppw)+2*nbdy;  % number of lateral elements
nZ = round(wZ/lambda*ppw)+2*nbdy;  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);

dY = c0/omega0*2*pi/ppw
dZ = c0/omega0*2*pi/ppw
dT = dY/c0*cfl;

%id = fopen(['c.dat'],'wb'); fwrite(fid,extendMap(cmap,nbdy),'float'); fclose(fid);
%fid = fopen(['A.dat'],'wb'); fwrite(fid,extendMap(Amap,nbdy),'float'); fclose(fid);
%fid = fopen(['rho.dat'],'wb'); fwrite(fid,extendMap(rhomap,nbdy),'float'); fclose(fid);
Nmap=(1+boveramap/2)./(rhomap.*cmap.^4);
%fid = fopen(['N.dat'],'wb'); fwrite(fid,extendMap(Nmap,nbdy),'float'); fclose(fid);

ncoords = size(incoords,1);
ncoordsout = size(outcoords,1);
incoords = incoords+nbdy;
outcoords = outcoords+nbdy;

try6_rebuild2(nY,nZ,nT,single(dY),single(dZ),single(dT),single(c0),ncoords,ncoordsout,int32(incoords),int32(outcoords),single(extendMap(cmap',nbdy)),single(extendMap(rhomap',nbdy)),single(extendMap(Amap',nbdy)),single(extendMap(Nmap',nbdy)),single(icmat'));
