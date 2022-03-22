%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% FIRST WRITTEN: 2018-06-21
% LAST MODIFIED: 2022-03-03
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % average speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wX = 2e-2;         % width of simulation field (m)
wY = 3e-2;         % depth of simulation field (m)
duration = wY*2.3/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;          % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition
ppp = ppw/cfl;     % points per period
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi; % wavelength (m)
nX = round(wX/lambda*ppw);  % number of lateral elements
nY = round(wY/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl); % number of time points
dX = c0/omega0*2*pi/ppw % step size in x
dY = dX; % step size in y (please keep step sizes the same)
dT = dX/c0*cfl; % step size in time

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
inmap(:,1:8)=ones(nX,8);
imagesc(inmap'), axis equal, axis tight
incoords = mapToCoords(inmap); % note zero indexing for compiled code
plot(incoords(:,1),incoords(:,2),'.')
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
fcen=[round(nX/2) round(nY/1.3)]; % center of focus
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat=repmat(icvec,size(incoords,1)/8,1);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' repmat(icvec,size(incoords,1)/8,1)']';
end
imagesc(icmat)
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
modX=1; modY=1;
outcoords = coordsMatrix(nX,nY,modX,modY);
plot(outcoords(:,1),outcoords(:,2),'.')
outcoords(:,3)=0; % label zero for total field
outcoords(find(outcoords(:,2)==9),3)=1; % label 1 for transducer surface
idc=find(outcoords(:,3)==0);
idxducer=find(outcoords(:,3)==1);
plot(outcoords(idc,1),outcoords(idc,2),'.'),hold on
plot(outcoords(idxducer,1),outcoords(idxducer,2),'r.'), hold off
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap = ones(nX,nY)*1000; % density map (kg/m^3)
Amap = ones(nX,nY)*0.0;    % attenuation map (dB/MHz/cm)
betamap = ones(nX,nY)*0.0;    % nonlinearity map 
%% Let's add a surface %%%%%%%%%%%%%%%%%%%%%%%%%%
cmap(:,round(nY/2)+8:end)=0.95*c0; % surface
imagesc(cmap')
prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end

pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))

px=pxducer(:,round(end/2)); % rf line at center of transducer
plot(px)
[val idt0]=max(abs(hilbert(px))); % idt0 is the time pixel at which the transducer emission envelope is maximal

% convert time into space
tyaxis=((1:size(pxducer,1))-idt0)*dT*c0/2;
plot(tyaxis,px), grid on
%compare to surface location
round(nY/2)*dY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's look at axial resolution by putting another surface right behind the first one%
% What do you think will happen if you change the number of cycles?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add a surface %%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
cmap(:,round(nY/2)+8:end)=0.95*c0; % first surface
cmap(:,round(nY/2)+8+ppw:end)=c0; % second surface
imagesc(cmap')
prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end

pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))

px=pxducer(:,round(end/2)); % rf line at center of transducer
plot(px)
[val idt0]=max(abs(hilbert(px))); % idt0 is the time pixel at which the transducer emission envelope is maximal

% convert time into space
tyaxis=((1:size(pxducer,1))-idt0)*dT*c0/2;
plot(tyaxis,px), grid on
%compare to surface location
round(nY/2)*dY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's add the scatterer in and let's focus to the scatterer%%
%% What would happen if you didn't focus on transmit on the scattererr? %%%
foc=round(nY/1.3);
fcen=[round(nX/2) foc]; % center of focus on axis
fcen=[round(nX/10) foc]; % center of focus off-axis

cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap=ones(nX,nY)*1000;
rhomap(fcen(1)-1:fcen(1)+1,fcen(2)-1:fcen(2)+1)=500; % scatterer
imagesc(rhomap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end

pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))

px=pxducer(:,round(end/2)); % rf line at center of transducer
plot(px)
[val idt0]=max(abs(hilbert(px))); % idt0 is the time pixel at which the transducer emission envelope is maximal

% convert time into space
tyaxis=((1:size(pxducer,1))-idt0)*dT*c0/2;
plot(tyaxis,px), grid on
%compare to scatterer location
round(nY/1.3)*dY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What happens if we move our scatterer off-axis?
%% let'c change fcen and find out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What happens if we use a plane wave transmit? 
foc=round(nY/1.3);
fcen=[round(nX/2) foc*1000]; % center of focus on axis

fcen_scat=[round(nX/2) foc]; % center of focus on axis
fcen_scat=[round(nX/10) foc]; % center of focus off-axis


cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap=ones(nX,nY)*1000;
rhomap(fcen_scat(1)-1:fcen_scat(1)+1,fcen_scat(2)-1:fcen_scat(2)+1)=500; % scatterer
imagesc(rhomap'), axis equal, axis tight

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing')
toc
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end

pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))

px=pxducer(:,round(end/2)); % rf line at center of transducer
plot(px)
[val idt0]=max(abs(hilbert(px))); % idt0 is the time pixel at which the transducer emission envelope is maximal

% convert time into space
tyaxis=((1:size(pxducer,1))-idt0)*dT*c0/2;
plot(tyaxis,px), grid on
%compare to scatterer location
round(nY/1.3)*dY

% different location on transducer
plot(tyaxis,px), hold on,
plot(tyaxis,pxducer(:,round(end/10))), hold off
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point targets are spherical reflectors
% Let's pick an aribitrary point in space and calculate the distance to the transducer surface

xducercoords = outcoords(idxducer,:);
i=nX/2; j=nY/2; % point in the middle
i=fcen_scat(1)+2*ppw; j=fcen_scat(2); % location of scatterer
i=fcen_scat(1); j=fcen_scat(2); % location of scatterer

subplot(2,1,1)
plot(i,j,'.'), hold on
plot(xducercoords(:,1),xducercoords(:,2),'r.'), hold off
axis equal, xlim([1 nX]), ylim([1 nY])

dd=zeros(size(xducercoords,1),1);
for n=1:size(xducercoords,1)
  subplot(2,1,1)
  plot(i,j,'.'), hold on
  plot(xducercoords(:,1),xducercoords(:,2),'r.')
  plot(xducercoords(n,1),xducercoords(n,2),'g.')
  hold off
  axis equal, xlim([1 nX]), ylim([1 nY])

  dd(n)=sqrt((xducercoords(n,1)-i).^2+(xducercoords(n,2)-j).^2);
  subplot(2,1,2)
  plot(dd)
  ylim([0 max(dd)])
  drawnow
end


%% Apply the delays to the RF data in three parts
% 1. idt0
% 2. range equation to calculate the arrival time
% 3. spherical delays for Huygen's principle

lat=0*(nY/2)*dX;
dep=(nY/2)*dY;

subplot(2,1,1), imagesc(powcompress(pxducer,1/3))
rftmp=pxducer;
rfrng=2*ppp;
rffoc=zeros(2*rfrng+1,size(xducercoords,1));
idt=idt0+(j-8)*dX/double(c0)/dT; % transmit offset + range equation
for n=1:size(xducercoords,1)
  m=round(idt+dd(n)/cfl-rfrng:idt+dd(n)/cfl+rfrng);
  rffoc(:,n)=pxducer(m,n); % the focused rf
  subplot(2,1,2), imagesc(rffoc), drawnow;
  rftmp(m,n)=NaN;
  subplot(2,1,1), imagesc(powcompress(rftmp,1/3)), drawnow;

end
    
%% Now do it for all points in the image you want to reconstruct
imagesc(powcompress(pxducer,1/3))
deps = 1e-3:lambda/8:nY*0.9*dX;
lats = -20e-3:lambda/8:20e-3;

xducercoords = outcoords(idxducer,:);
bm=zeros(length(lats),length(deps));
idps=cell(length(lats),length(deps));

fnumber=1;

idps=cell(length(lats),length(deps));
for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fcen=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
    idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
    dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
    idt=idt0+round(2*dep/double(c0)/(dT));
    idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
    idps{ii,jj}=idp;
  end
end

  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj)=sum(pxducer(idps{ii,jj}));
    end
  end

figure(1), clf
imagesc(lats*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm')))),[-40 0])
colormap gray, cbar=colorbar; title(cbar,'dB')
xlabel('mm'), ylabel('mm')
axis equal, axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Focus on transmit and focus on receive%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Walking aperture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx.nTx=40; % number of Tx events
tx.dep=nY/1.3*dY; % focal depth (m)
tx.bmw=lambda/2; % beamwidth (m), is this an integer multiple of dX? 

fcen=[round(nX/2) tx.dep/dX]; % center of focus on axis

t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
[icmat dd] = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/8,:),icvec,cfl);
for k=2:8
  t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
  icmat=[icmat' focusCoords(fcen(1),fcen(2),incoords((k-1)*size(incoords,1)/8+1:(k)*size(incoords,1)/8,:),icvec,cfl)']';
end
imagesc(icmat)

% instead of moving the transducer, we will move the maps
nXextend=nX+ceil(tx.bmw/dX*tx.nTx);
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmapextend = ones(nXextend,nY)*1540;   % speed of sound map (m/s)
rhomapextend = ones(nXextend,nY)*1000; % density map (kg/m^3)
Amapextend = ones(nXextend,nY)*0.0;    % attenuation map (dB/MHz/cm)
betamapextend = ones(nXextend,nY)*0.0;    % nonlinearity map 

scat_density=0.0015;
scats=rand(nXextend,nY);
scats(find(scats>scat_density))=0;
scats=scats/max(max(scats));
mean(mean(scats))
scats(:,1:10)=0; % don't put scatters inside your transducer
rhosr=0.0375; % scatterer impedance contrast 

rhomapextend=rhomapextend-scats*1000*rhosr;
imagesc(rhomapextend'), colorbar

for n=1:tx.nTx
  outdir=['/kulm/scratch/bmm890/txrx_' num2str(n) '/']
  eval(['!mkdir -p ' outdir]);

  cmap=cmapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
  rhomap=rhomapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
  Amap=Amapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
  betamap=betamapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
  
  imagesc(rhomap'), drawnow

  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  
  cwd=pwd; addpath(cwd);
  cd(outdir)
  prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
  eval('!./fullwave2_try6_nln_relaxing & ')
  cd(cwd);

end


%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 2e-3:lambda/8:nY*dY/1.1;
lats = 0;
xducercoords = outcoords(idxducer,:);
bm=zeros(length(lats),length(deps),tx.nTx);
idps=cell(length(lats),length(deps));
fnumber=1;

n=round(tx.nTx/2);
 outdir=['/kulm/scratch/bmm890/txrx_' num2str(n) '/']
ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))
px=pxducer(:,round(size(pxducer,2)/2));
[val idt0]=max(abs(hilbert(px)))

for n=1:tx.nTx
  outdir=['/kulm/scratch/bmm890/txrx_' num2str(n) '/']

  ncoordsout=size(outcoords,1);
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  while(nRun<nT-1)
    pause(0.1)
    nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  end
 pxducer = readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
 imagesc(powcompress(pxducer,1/3)), drawnow
 
  if(n==1)
    idps=cell(length(lats),length(deps));
    for ii=1:length(lats)
      lat=lats(ii);
      for jj=1:length(deps)
        dep=deps(jj);
        fcen=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
        idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
        dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
        idt=idt0+round(2*dep/double(c0)/(dT));
        idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
	idp=idp(find(idp<=size(pxducer,1)*size(pxducer,2)));
        idps{ii,jj}=idp;
      end
    end
  end

  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(pxducer(idps{ii,jj}));
    end
  end
end


%% PLOT THE BMODE IMAGE %%
figure(1)
n=1:tx.nTx; bws=((n-(tx.nTx+1)/2)*tx.bmw);
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight

% think about speed of sound, idt0, colorbar, sampling, ...

addpath /celerina/gfp/mfs/dumbmat/
figure(1)
n=1:tx.nTx; bws=((n-(tx.nTx+1)/2)*tx.bmw);
img=dbzero(abs(hilbert(squeeze(bm))));
img=interp2easy(img,4,1);
imagesc(bws*1e3,deps*1e3,img,[-40 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare to plane wave %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=round(tx.nTx/2);

outdir=['/kulm/scratch/bmm890/txrx_' num2str(n) '/']
eval(['!mkdir -p ' outdir]);

cmap=cmapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
rhomap=rhomapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
Amap=Amapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);
betamap=betamapextend(1+(n-1)*round(tx.bmw/dX):(n-1)*round(tx.bmw/dX)+nX,:);

imagesc(rhomap'), drawnow

eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);

cwd=pwd; addpath(cwd);
cd(outdir)
prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
eval('!./fullwave2_try6_nln_relaxing & ')
cd(cwd);

outdir=['/kulm/scratch/bmm890/txrx_' num2str(n) '/']
ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
pxducer = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))
px=pxducer(:,round(size(pxducer,2)/2));
[val idt0]=max(abs(hilbert(px)))

imagesc(powcompress(pxducer,1/3))
lats = bws;

bm=zeros(length(lats),length(deps));
idps=cell(length(lats),length(deps));

fnumber=1;

idps=cell(length(lats),length(deps));
for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fcen=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
    idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
    dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
    idt=idt0+round(2*dep/double(c0)/(dT));
    idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
    idps{ii,jj}=idp;
  end
end

  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj)=sum(pxducer(idps{ii,jj}));
    end
  end

figure(1), clf
imagesc(lats*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm')))),[-40 0])
colormap gray, cbar=colorbar; title(cbar,'dB')
xlabel('mm'), ylabel('mm')
axis equal, axis tight
