%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% FIRST WRITTEN: 2018-06-21
% LAST MODIFIED: 2022-03-03
% Launch Fullwave 2 code, easy matlab wrapper
% Easy plane wave beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % average speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wX = 3e-2;         % width of simulation field (m)
wY = 5e-2;         % depth of simulation field (m)
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
nangles=11; 
dtheta=1.75*pi/180;
rho0=1000;

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
for m=1:1:floor(nY*dY*1e2)
  rhomap(round(nX/2),round(m*1e-2/dY))=0.65*rho0; % point
end
imagesc(rhomap')

if(0) % speckle
  scat_density=0.015;
  scats=rand(nX,nY);
  scats(find(scats>scat_density))=0;
  scats=scats/max(max(scats));
  mean(mean(scats))
  scats(:,1:10)=0; % don't put scatterers inside your transducer
  rhosr=0.0375; % scatterer impedance contrast 
  rhomap=rhomap-scats*rho0*rhosr;
end


%[m.c m.rho m.A m.N] = img2fieldFlatten('r102gh.tif',dX,dY,c0,rho0);

c=cmap; rho=rhomap; A=Amap; beta=betamap;
for n=1:nangles
    theta=(n-(nangles+1)/2)*dtheta

    fcen=[round(1e6/dY)*sin(theta) round(1e6/dY)]; % center of focus
    t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
    icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
    plot(icvec)
    icmat=repmat(icvec,size(incoords,1)/8,1);
    for k=2:8
      t=t-dX/c0; icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
      icmat=[icmat' repmat(icvec,size(incoords,1)/8,1)']';
    end
    icmat = focusCoords(fcen(1),fcen(2),incoords,icvec,cfl);
    imagesc(icmat), drawnow
  
    cwd=pwd; addpath(cwd);
    outdir=['/kulm/scratch/lesion' num2str(n)]; eval(['!mkdir -p ' outdir]); 
    eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
    cd(outdir)
 
    prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)
  if(n<nangles)
    eval('!./fullwave2_try6_nln_relaxing &')
  elseif(n==nangles)
    eval('!./fullwave2_try6_nln_relaxing')
  end
  cd(cwd);
end


%% propagation movie %%
n=round(nangles/2);
outdir=['/kulm/scratch/lesion' num2str(n) '/']
idc=find(outcoords(:,3)==0);

ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout
genout = readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idc);
p = reshape(genout,size(genout,1),nY2-1,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

%% use powcompress to view small amplitude fields %%
for i=1:20:size(p,1)
  imagesc(powcompress(squeeze(p(i,:,:)),1/3)), title(num2str(i)), drawnow
end



%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 1e-3:lambda/8:nY*dY;
lats = -wX/2:lambda/8:wX/2;

idxducer=find(outcoords(:,3)==1);
xducercoords = outcoords(idxducer,:);
ncoordsout=size(outcoords,1);

n=round(nangles/2);
outdir=['/kulm/scratch/lesion' num2str(n) '/']
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;

pxducer = readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
imagesc(powcompress(pxducer,1/3))
px=pxducer(:,round(size(pxducer,2)/2));
[val idt0]=max(abs(hilbert(px)))

bm=zeros(length(lats),length(deps),nangles);
idps=cell(length(lats),length(deps));
fnumber=1;

for n=1:nangles
    theta=(n-(nangles+1)/2)*dtheta
     
  outdir=['/kulm/scratch/lesion' num2str(n) '/']
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;

  pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idxducer);
  figure(3), imagesc(powcompress(pxducer,1/4)), drawnow
  px=pxducer(:,round(size(pxducer,2)/2));
  [val idt0]=max(abs(hilbert(px)))

    idps=cell(length(lats),length(deps));
    for ii=1:length(lats)
      lat=lats(ii);
      for jj=1:length(deps)
        dep=deps(jj);
        fcen=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
        idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
        dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
        idt=idt0+round(2*dep/double(c0)/(dT));
        %idt=idt0+round(dep/double(c0)/(dT)+(dep*cos(theta)+lat*sin(theta))/double(c0)/dT);
         idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
	idp=idp(find(idp>0 & idp<=size(pxducer,1)*size(pxducer,2)));
        idps{ii,jj}=idp;
      end
    end

  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(pxducer(idps{ii,jj}));
    end
  end

  figure(1)
  imagesc(lats,deps,dbzero(abs(hilbert(bm(:,:,n)'))),[-50 0]);
   colormap gray, cbar=colorbar; title(cbar,'dB')
  xlabel('mm'), ylabel('mm')
  axis equal, axis tight

  
  figure(2)
  imagesc(lats,deps,dbzero(abs(hilbert(mean(bm(:,:,1:n),3)'))),[-50 0]);
  colormap gray, cbar=colorbar; title(cbar,'dB')
  xlabel('mm'), ylabel('mm')
  axis equal, axis tight

end

%% PLOT THE BMODE IMAGE %%

for n=1:nangles
  figure(2)
  imagesc(lats,deps,dbzero(abs(hilbert(mean(bm(:,:,1:n),3)'))),[-50 0]);
  colormap gray, cbar=colorbar; title(cbar,'dB')
  xlabel('mm'), ylabel('mm')
  axis equal, axis tight
  drawnow
end





