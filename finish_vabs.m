lambda = c0/omega0*2*pi;
dY = lambda/ppw;
dT = dY/c0*cfl;
dZ = dY;

tMin = -8*ncycles*pi/(omega0);
N = -beta/(rho0*c0^4)

nZ = round(dep/dZ);
dT = dZ/c0/5;
nT = round(nZ*dZ/c0/dT*2.5);

nYducer = round(ay/dY);
nY = round(wY/dY);

%modY = round(0.5e-4/dY);
%modZ = round(0.5e-4/dZ);
%modT = round(1/(100e6*dT));

nY*nZ*nT*4/1e9/modY/modZ/modT

nY2 = ceil(nY/modY); dY2 = dY*modY;
nZ2 = ceil(nZ/modZ); dZ2 = dZ*modZ;
nT2 = ceil(nT/modT)-1; dT2 = dT*modT;

zaxis = (0:nZ-1)*dZ;
yaxis = (0:nY-1)*dY; 
yaxis = yaxis-mean(yaxis); 

zaxis2 = (0:nZ2-1)*dZ2;
yaxis2 = (0:nY2-1)*dY2; 
yaxis2 = yaxis2-mean(yaxis2); 
