%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: APR 17, 2017
% Fullwave imaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wY = 2.4e-2;         % width of simulation field (m)
wZ = 6e-2;         % depth of simulation field (m)
p0 = 1e5; % pressure in Pa
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;           % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nY = round(wY/lambda*ppw);  % number of lateral elements
nZ = round(wZ/lambda*ppw);  % number of depth elements
dY = c0/omega0*2*pi/ppw
dZ = c0/omega0*2*pi/ppw
duration = nZ*dZ/c0*2.5;  % duration of simulation (s)
nT = round(duration*c0/lambda*ppw/cfl);
dT = dY/c0*cfl;
foc=round(nZ/1.3);
fcen=[round(nY/2) foc]; % center of focus
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nY,nZ); 
inmap(:,1) = ones(nY,1); inmap(:,2) = ones(nY,1); inmap(:,3) = ones(nY,1);
incoords = mapToCoords(inmap);
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
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap = zeros(nY,nZ);  outmap(:,3) = ones(nY,1);
outcoords = mapToCoords(outmap);
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0=1000;
A0 = 3/(20/log(10));
N0 = -tissue.beta/(rho0*c0^4);
fnumber=foc/nY;
beamwidth=round(lambda*fnumber/dY);
nlines=11;
[m.c m.rho m.A m.N] = img2fieldFlatten('r102gh.tif',dY,dZ,c0,rho0);
nYextend=size(m.c,1);
res_cell = rescell2d(c0,omega0,foc*dY,wY,ncycles,dY,dZ);
num_scat = 20;
scat_size = dZ;
scat_size_samp = round(scat_size/dY);
cscat_lesion = generate_c_scat(1,1,num_scat/res_cell, scat_size_samp, nYextend, nZ)-1;
idl=circleIdx(size(cscat_lesion),[nYextend/2 foc],3e-3/dY);
cscat_lesion(idl)=0;
imagesc(cscat_lesion'), colorbar
csr=0.05; % scatterer impedance contrast 
for n=1:nlines
  orig = [round(1.5e-2/dY-(n-(nlines+1)/2)*beamwidth/4) 1]
  c = chopField(m.c,c0,orig,nY,nZ);
  rho = chopField(m.rho,rho0,orig,nY,nZ);
  A = chopField(m.A,A0,orig,nY,nZ);
  N = chopField(m.N,N0,orig,nY,nZ); 
  cs = cscat_lesion(round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nY/2):round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nY/2)+nY-1,:);
  bovera=N*0-2;

  c=c*0+c0;
  rho=rho*0+rho0;
  A=A*0+A0;
  
  c=c-cs*c0*csr;
  imagesc(c'), drawnow

  cwd=pwd; addpath(cwd);
  outdir=['lesion' num2str(n)]; eval(['!mkdir -p ' outdir]); 
  eval(['!cp try6_nomex ' outdir]);
  cd(outdir)
  launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,c',rho',A',bovera',incoords,outcoords,icmat);
  if(n<nlines)
    eval('!./try6_nomex &')
  elseif(n==nlines)
    eval('!./try6_nomex ')
  end
  cd(cwd);
end
%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 10e-3:0.125e-3/4:foc*dZ*1.1;
lats = 0;
xducercoords = outcoords;
bm=zeros(length(lats),length(deps),nlines);
idps=cell(length(lats),length(deps));

n=round(nlines/2);
outdir=['lesion' num2str(n) '/']
ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1));
imagesc(powcompress(pxducer,1/4))
px=pxducer(:,round(size(pxducer,2)/2));
[val idt0]=max(abs(hilbert(px)))

for n=1:nlines
  outdir=['lesion' num2str(n) '/']
  ncoordsout=size(outcoords,1);
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  while(nRun~=nT)
    pause(0.1)
    nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  end
  pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1));
  imagesc(powcompress(pxducer,1/4))
  px=pxducer(:,round(size(pxducer,2)/2));
  [val idt0]=max(abs(hilbert(px)))

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
n=1:nlines; bws=((n-(nlines+1)/2)*beamwidth/4)*dY;
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray
xlabel('mm'), ylabel('mm')
axis equal, axis tight

figure(2)
n=6;
cs = cscat_lesion(round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nY/2):round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nY/2)+nY-1,:);

imagesc(((1:nY)-nY/2)*dY*1e3,(1:nZ)*dY*1e3,cs')
axis equal, axis tight
axis([bws(1) bws(end) deps(1) deps(end)]*1e3)
xlabel('mm'), ylabel('mm')
