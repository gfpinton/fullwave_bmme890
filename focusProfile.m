function [dd] = focusProfile (fcen,coords,cfl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% focal profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dd=zeros(size(coords,1),1);
for i=1:length(fcen)
  dd = dd + (coords(:,i)-fcen(i)).^2;
end
dd=sqrt(dd);
%  dd = sqrt((coords(:,1)-fcen(1)).^2+(coords(:,2)-fcen(2)).^2+(coords(:,3)-fcen(3)).^2);
dd = round(dd/cfl);
dd = dd-min(dd);

