function [foc_prof] = focprof(xaxis, yaxis, dT, c0, dx, dy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% establish a focusing profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X Y] = meshgrid(xaxis,yaxis);
foc_prof = sqrt(X.^2+dx^2)/c0-dx/c0+sqrt(Y.^2+dy^2)/c0-dy/c0;
% convert to samples
foc_prof = foc_prof/dT;
foc_prof = round(foc_prof);
%  figure(1),imagesc(foc_prof), colorbar
% pxducer = zeros(nY,nX,500);

