%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2018-06-21
% LAST MODIFIED: 2022-04-12
% Launch Fullwave 2 code, easy matlab wrapper
% Time reversal focusing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % average speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wX = 4e-2;         % width of simulation field (m)
wY = 3e-2;         % depth of simulation field (m)
duration = wY*2.3/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;          % number of points per spatial wavelength
cfl = 0.3;         % Courant-Friedrichs-Levi condition
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
foc=round(nY/1.3); rad=10% 
fcen=[round(nX/2) foc]; % center of focus
inmap = zeros(nX,nY); 
for i=1:nX
    for j=1:nY
        if(sqrt((i-fcen(1))^2+(j-fcen(2))^2)<=rad)
            inmap(i,j)=1;
        end
    end
end
incoords = mapToCoords(inmap);
plot(incoords(:,1),incoords(:,2),'.')
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
fcen=[round(nX/2) round(nY/1.3)]; % center of focus
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat=repmat(icvec,size(incoords,1),1);
imagesc(icmat)
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
modX=1; modY=1;
outcoords = coordsMatrix(nX,nY,modX,modY);
outcoords(:,3)=1; % Full matrix output label
idc=find(outcoords(:,2)<4); outcoords(idc,3)=2; % transducer surface label
plot(outcoords(:,1),outcoords(:,2),'.'),hold on
plot(outcoords(idc,1),outcoords(idc,2),'r.'), hold off
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nX,nY)*1540;   % speed of sound map (m/s)
rhomap = ones(nX,nY)*1000; % density map (kg/m^3)
Amap = ones(nX,nY)*0.0;    % attenuation map (dB/MHz/cm)
betamap = ones(nX,nY)*0.0;    % nonlinearity map
if(0) % simple distributed heterogeneities
cmap(round(nX/4):round(nX/2),round(nY/3):round(nY/2))=c0*0.8;
cmap(round(nX/2):round(nX),round(nY/3):round(nY/2))=c0*0.9;
cmap(round(nX/10):round(nX/2),round(nY/2):round(3*nY/4))=c0*0.85;
end

if(0) % complex distributed heterogeneites
for n=1:100
  r=rand(1)*min([nX nY])/10
  idl=circleIdx(size(cmap),[nX*rand(1) nY/2*rand(1)],r);
  cmap(idl)=c0*(rand(1)-0.5)/3+c0;
  imagesc(cmap'),colorbar, drawnow
end
end
if(0) % human  skull
load skullmap2d
cskull=interp2easy(skullmap2d.c,skullmap2d.dX/dX,skullmap2d.dX/dX);
cskull=cskull(round(109.2e-3/dX)-round(nX/2)+1:round(109.2e-3/dX)-round(nX/2)+nX,round(25.7e-3/dX):round(25.7e-3/dX)+nY-1);
cmap=cskull;
end

imagesc(cmap'), axis equal, axis tight
%%% write files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,outcoords,icmat)
tic
eval('!./fullwave2_try6_nln_relaxing &')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)
nX2=length(1:modX:nX); nY2=length(1:modY:nY);
pause(1)
nRun=sizeOfFile('genout.dat')/4/ncoordsout
while nRun<nT-1
  genout = readGenoutSlice(['genout.dat'],nRun-1,size(outcoords,1));
  pslice = reshape(genout,size(genout,1),nY2,nX2);
  imagesc(squeeze(pslice(end,:,:))), colorbar, drawnow
  nRun=sizeOfFile('genout.dat')/4/ncoordsout
end

nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:20:size(p,1)
  imagesc(squeeze(p(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW TIME REVERSE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idc=find(outcoords(:,3)==2);
pxducer1=genout(:,idc);
imagesc(powcompress(pxducer1,1/4)) % transducer plotted on a compressed scale

px1=pxducer1(:,round(size(pxducer1,2)/2));
[val idt0]=max(abs(hilbert(px1)))
plot(px1)
icmat_tr=pxducer1(idt0-10*ppw:idt0+40*ppw,:)';
icmat_tr=flipdim(icmat_tr,2); % time reversal
icmat_tr=icmat_tr/(max(max(icmat_tr)))*p0;
imagesc(icmat_tr);
icmat_tr(:,end+1:nT)=0;

incoords_tr = outcoords(idc,:);
%%% write files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep_fullwave2_try6_nln_relaxing9(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords_tr,outcoords,icmat_tr)
tic
eval('!./fullwave2_try6_nln_relaxing &')
toc
pause(1)
nRun=sizeOfFile('genout.dat')/4/ncoordsout
while nRun<nT-1
  genout = readGenoutSlice(['genout.dat'],nRun-1,size(outcoords,1));
  pslice = reshape(genout,size(genout,1),nY2,nX2);
  imagesc(squeeze(pslice(end,:,:))), colorbar, drawnow
  nRun=sizeOfFile('genout.dat')/4/ncoordsout
end

nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),nY2,nX2);
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:20:size(p,1)
  imagesc(squeeze(p(i,:,:)), [-1 1]*p0), title(num2str(i)), drawnow
end
imagesc(squeeze(sum(p.^2)))
imagesc(dbzero(squeeze(sum(p.^2)))/2), caxis([-40 0]), colorbar
