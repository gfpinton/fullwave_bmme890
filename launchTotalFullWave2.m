function [dY dZ] = launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap,rhomap,Amap,boveramap,incoords,outcoords,icmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Write simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbdy = 40; % number of boundary points for PML

lambda = c0/omega0*2*pi;
lambda = c0/omega0*2*pi;
nY = round(wY/lambda*ppw)+2*nbdy;  % number of lateral elements
nZ = round(wZ/lambda*ppw)+2*nbdy;  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);

dY = c0/omega0*2*pi/ppw
dZ = c0/omega0*2*pi/ppw
dT = dY/c0*cfl;

nTic=size(icmat,2);

fid = fopen(['c.dat'],'wb'); fwrite(fid,extendMap(cmap,nbdy),'float'); fclose(fid);
fid = fopen(['A.dat'],'wb'); fwrite(fid,extendMap(Amap,nbdy),'float'); fclose(fid);
fid = fopen(['rho.dat'],'wb'); fwrite(fid,extendMap(rhomap,nbdy),'float'); fclose(fid);
fid = fopen(['N.dat'],'wb'); fwrite(fid,extendMap((1+boveramap/2)./(rhomap.*cmap.^4),nbdy),'float'); fclose(fid);

ncoords = size(incoords,1);
ncoordsout = size(outcoords,1);
incoords = incoords+nbdy;
outcoords = outcoords+nbdy;

writeCoords('icc.dat',incoords);
writeCoords('outc.dat',outcoords);
writeIC('icgen.dat',icmat);

writeVabs('float',dY,'dY',dZ,'dZ',dT,'dT',c0,'c0');
writeVabs('int',nY,'nY',nZ,'nZ',nT,'nT',ncoords,'ncoords',ncoordsout,'ncoordsout',nTic,'nTic');
%outdir = './';
%writeVabs('char',outdir,'outdir');



%try3(nY,nZ,nT,dY,dZ,dT,c0,ncoords,ncoordsout);
