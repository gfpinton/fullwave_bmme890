%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Launch Fullwave code, easy matlab wrapper
% Beamforming and point spread function (PSF) calculations
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

%foc=round(nZ/1.3);
%fcen=[round(nY/2) foc]; % center of focus
%z1=0.95;
%z2=0.925;
%cmap(round(nY/2)-1:round(nY/2)+1,foc-1)=z1*c0;
%cmap(round(nY/2)-1:round(nY/2)+1,foc)=z2*c0;
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nY,nZ); 
inmap(round(nY/4),round(nZ/3))=1;
inmap(round(nY*3/4),round(nZ*2/3))=1;
imagesc(inmap')
incoords = mapToCoords(inmap);
%%% Generate initial conditions based on input coordinates %%%%%%
icmat=zeros(size(incoords,1),nT);

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec1 = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec1), hold all
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi -round(nZ/3)*dZ/c0;
icvec2 = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;


icmat=[icvec2; icvec1;];
imagesc(icmat)



%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap=zeros(nY,nZ); modY=1; modZ=1;
[modidy modidz]=meshgrid(1:modY:nY,1:modZ:nZ);
outmap(modidy,modidz)=1;
imagesc(outmap');
outcoords=mapToCoords(outmap);
outcoords(:,3)=1;
outcoords=[outcoords' [incoords 2*ones(size(incoords,1),1)]']';
idc=find(outcoords(:,3)==1);
plot(outcoords(idc,1),outcoords(idc,2),'.')
hold on
idc=find(outcoords(:,3)==2);
plot(outcoords(idc,1),outcoords(idc,2),'r.')
hold off
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
for i=1:1:size(p1,1)
  imagesc(squeeze(p1(i,:,:))', [-1 1]*p0), title(num2str(i)), drawnow
end

pxducer=p1(:,:,3);



idc=find(outcoords(:,3)==2);
pxducer1=genout(:,idc);
pxducer1=pxducer1(:,round(size(pxducer1,2)/3*2+1:size(pxducer1,2)));
imagesc(powcompress(pxducer1,1/4)) % transducer plotted on a compressed scale
px1=pxducer1(:,round(size(pxducer1,2)/2));
[val idt0]=max(abs(hilbert(px1)))
plot(px1)
idt=idt0+round(2*foc*dZ/c0/dT);
plot(px1(idt-3*round(ppw/cfl):idt+2*round(ppw/cfl))) % note slightly offcenter due to dispersion

%%% GENERATE PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 10e-3:0.125e-3/4:foc*dZ*1.1;
lats = -5e-3:0.25e-3/4:5e-3;
bm1=zeros(length(lats),length(deps));
fnumber=1;
idc=find(outcoords(:,3)==2);idc=idc(round(length(idc)/3*2+1:length(idc)));
xducercoords = outcoords(idc,:);

dX=dY; 
pxducer=pxducer1;
idps=cell(length(lats),length(deps));
for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fcen=round([lat/dX+mean(xducercoords(:,1)) dep/dX ]);
    idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber); % constant Fnumber
    dd=focusProfile(fcen,xducercoords(idx,:),dT/dX*c0);
    idt=idt0+round(2*dep/double(c0)/(dT));
    idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
    idps{ii,jj}=idp;
  end
end
bm1=zeros(length(lats),length(deps),'single');
tic
for ii=1:length(lats)
  for jj=1:length(deps)
    bm1(ii,jj)=sum(pxducer1(idps{ii,jj}));
  end
end
toc
imagesc(lats*1e3,deps*1e3,dbzero(abs(hilbert(bm1'))),[-45 0]), colormap gray, drawnow
xlabel('mm'),ylabel('mm')
axis equal, axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD ABDOMINAL MAP %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0=1000;
A0 = 3/(20/log(10));
N0 = -tissue.beta/(rho0*c0^4);
[m.c m.rho m.A m.N] = img2fieldFlatten('r102gh.tif',dY,dZ,c0,rho0);
imagesc(m.c'), colormap jet, axis equal, axis tight % total abdomen
size(m.c,1)*dY
orig = [round(1.5e-2/dY) 1];
c = chopField(m.c,c0,orig,nY,nZ);
imagesc(c') % reduce abdominal map to simulation size
rho = chopField(m.rho,rho0,orig,nY,nZ);
A = chopField(m.A,A0,orig,nY,nZ);
N = chopField(m.N,N0,orig,nY,nZ); 
bovera=N*0-2;
c(round(nY/2)-1:round(nY/2)+1,foc-1)=z1*c0;
c(round(nY/2)-1:round(nY/2)+1,foc)=z2*c0;
imagesc(c') % add scatterer
%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,c',rho',A',bovera',incoords,outcoords,icmat);
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
for i=1:5:size(p2,1)
  imagesc(squeeze(p2(i,:,:))', [-1 1]*p0), title(num2str(i)), drawnow
end

cmod=c(1:modY:end,1:modZ:end); cmod=cmod-min(min(cmod)); cmod=cmod/max(max(cmod));
for i=1:5:size(p2,1)
  imagesc(powcompress(squeeze(p2(i,:,:))'+cmod'*p0/20,1/3), powcompress([-1 1]*p0,1/3)), title(num2str(i)), drawnow
end

idc=find(outcoords(:,3)==2);
pxducer2=genout(:,idc);
pxducer2=pxducer2(:,round(size(pxducer2,2)/3*2+1:size(pxducer2,2)));
imagesc(powcompress(pxducer2,1/3))

bm2=zeros(length(lats),length(deps),'single');
tic
for ii=1:length(lats)
  for jj=1:length(deps)
    bm2(ii,jj)=sum(pxducer2(idps{ii,jj}));
  end
end
toc

figure(2), clf
imagesc(lats*1e3,deps*1e3,dbzero(abs(hilbert(bm2'))),[-45 0]), colormap gray, drawnow
xlabel('mm'),ylabel('mm')
axis equal, axis tight
cbar=colorbar; title(cbar,'dB')

figure(1), clf
imagesc(lats*1e3,deps*1e3,dbzero(abs(hilbert(bm1'))),[-45 0]), colormap gray, drawnow
xlabel('mm'),ylabel('mm')
axis equal, axis tight
cbar=colorbar; title(cbar,'dB')
