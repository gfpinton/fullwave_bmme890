%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: APRIL 13, 2018
% LAST MODIFIED: APRIL 13, 2018
% Launch Fullwave code, easy matlab wrapper
% Time reversal focusing
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wY = 2e-2;         % width of simulation field (m)
wZ = 3e-2;         % depth of simulation field (m)
duration = 2.0*wZ/c0;  % duration of simulation (s)
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

cmap(round(nY/4):round(nY/2),round(nZ/3):round(nZ/2))=c0*0.8;
cmap(round(nY/2):round(nY),round(nZ/3):round(nZ/2))=c0*0.9;
cmap(round(nY/10):round(nY/2),round(nZ/2):round(3*nZ/4))=c0*0.85;
%cmap(round(nY/3):round(nY/2),round(nZ/3):round(nZ/2))=c0*0.8;
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
foc=round(nZ/1.3); rad=10% 
fcen=[round(nY/2) foc]; % center of focus
inmap = zeros(nY,nZ); 
%inmap(:,1) = ones(nY,1); inmap(:,2) = ones(nY,1); inmap(:,3) = ones(nY,1);
for j=1:nY
    for k=1:nZ
        if(sqrt((j-fcen(1))^2+(k-fcen(2))^2)<=rad)
            inmap(j,k)=1;
        end
    end
end

incoords = mapToCoords(inmap);
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold all
icmat=ones(size(incoords,1),1)*icvec;
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
for i=1:10:size(p1,1)
      imagesc(squeeze(p1(i,:,:))', [-1 1]*p0), title(num2str(i)), drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW TIME REVERSE %%
idc=find(outcoords(:,2)<=2);
pxducer1=genout(:,idc);
%pxducer1=pxducer1(:,round(size(pxducer1,2)/3*2+1:size(pxducer1,2)));
imagesc(powcompress(pxducer1,1/4)) % transducer plotted on a compressed scale
px1=pxducer1(:,round(size(pxducer1,2)/2));
[val idt0]=max(abs(hilbert(px1)))
plot(px1)
%idt=idt0+round(2*foc*dZ/c0/dT);
%plot(px1(idt-3*round(ppw/cfl):idt+2*round(ppw/cfl))) % note slightly offcenter due to dispersion

icmat_tr=pxducer1(idt0-10*ppw:idt0+20*ppw,:)';
icmat_tr=flipdim(icmat_tr,2); % time reversal
icmat_tr=icmat_tr/(max(max(icmat_tr)))*p0;
imagesc(icmat_tr);
icmat_tr(:,end+1:nT)=0;

incoords_tr = outcoords(idc,:);

%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords_tr,outcoords,icmat_tr);
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
for i=1:5:size(p1,1)
  imagesc(squeeze(p1(i,:,:))', [-1 1]*p0), title(num2str(i)), drawnow
end
