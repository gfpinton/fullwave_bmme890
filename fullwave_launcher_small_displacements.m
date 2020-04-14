%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2016
% LAST MODIFIED: April 11, 2017
% Launch Fullwave code, easy matlab wrapper
% Small displacements
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wY = 2e-2;         % width of simulation field (m)
wZ = 3e-2;         % depth of simulation field (m)
duration = 40e-6;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 15;           % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nY = round(wY/lambda*ppw);  % number of lateral elements
nZ = round(wZ/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dY = c0/omega0*2*pi/ppw
dZ = c0/omega0*2*pi/ppw
dT = dY/c0*cfl;
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nY,nZ)*1540;   % speed of sound map (m/s)
rhomap = ones(nY,nZ)*1000; % density map (kg/m^3)
Amap = ones(nY,nZ)*0.0;    % attenuation map (dB/MHz/cm)
boveramap = -2*ones(nY,nZ);    % nonlinearity map 

foc=round(nZ/1.3);
fcen=[round(nY/2) foc]; % center of focus
z1=0.95;
z2=0.945;
cmap(round(nY/2)-1:round(nY/2)+1,foc-1)=z1*c0;
cmap(round(nY/2)-1:round(nY/2)+1,foc)=z2*c0;
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nY,nZ); 
inmap(:,1) = ones(nY,1); inmap(:,2) = ones(nY,1); inmap(:,3) = ones(nY,1);
incoords = mapToCoords(inmap);
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap=zeros(nY,nZ);
[modidy modidz]=meshgrid(1:2:nY,1:2:nZ);
outmap(modidy,modidz)=1;
outcoords=mapToCoords(outmap);
outcoords(:,3)=1;
outcoords=[outcoords' [incoords 2*ones(size(incoords,1),1)]']';
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold all
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold off
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);
imagesc(icmat)
%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords,outcoords,icmat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
!./try6_nomex
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
idc=find(outcoords(:,3)==1);
p1 = reshape(genout(:,idc),size(genout,1),size(modidy,2),size(modidy,1));
imagesc(squeeze(p1(end,:,:))), colorbar
%for i=1:size(p1,1)
%  imagesc(squeeze(p1(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
%end

idc=find(outcoords(:,3)==2);
pxducer1=genout(:,idc);
pxducer1=pxducer1(:,round(size(pxducer1,2)/3*2+1:size(pxducer1,2)));
imagesc(powcompress(pxducer1,1/4))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFY CMAP AND RELAUNCH %%
cmap = ones(nY,nZ)*1540;   % speed of sound map (m/s)
cmap(round(nY/2)-1:round(nY/2)+1,foc-1)=z2*c0;
cmap(round(nY/2)-1:round(nY/2)+1,foc)=z1*c0;
%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords,outcoords,icmat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
!./try6_nomex
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
idc=find(outcoords(:,3)==1);
p2 = reshape(genout(:,idc),size(genout,1),size(modidy,2),size(modidy,1));
imagesc(squeeze(p2(end,:,:))), colorbar

idc=find(outcoords(:,3)==2);
pxducer2=genout(:,idc);
pxducer2=pxducer2(:,round(size(pxducer2,2)/3*2+1:size(pxducer2,2)));
imagesc(powcompress(pxducer2,1/4))


plot(pxducer1(:,round(size(pxducer1,2)/2))-pxducer2(:,round(size(pxducer2,2)/2)))

px1=pxducer1(:,round(size(pxducer1,2)/2));
px2=pxducer2(:,round(size(pxducer2,2)/2));

[val idt0]=max(abs(hilbert(px1)))
idt=idt0+round(2*foc*dZ/c0/dT);
plot(px1(idt-round(3*ppw/cfl):idt+round(2*ppw/cfl)))
hold on
plot(px2(idt-round(3*ppw/cfl):idt+round(2*ppw/cfl)),'r')
hold off
grid on

