function [c] = anechoic_lesion(c,c0,nY,nZ,dY,dZ,pos,radius)
%pos = [0 dx]; radius = [2.5e-3 2.5e-3];

yaxis = (0:nY-1)*dY;
yaxis = yaxis-mean(yaxis);
zaxis = (0:nZ-1)*dZ;
for j=1:nY
  for k=1:nZ
    if(((yaxis(j)-pos(1))/radius(1))^2 + ((zaxis(k)-pos(2))/radius(2))^2 < 1)
      c(j,k) = c0;
    end
  end
end
