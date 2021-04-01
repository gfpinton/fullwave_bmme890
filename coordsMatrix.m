function [coords] = coordsMatrix(nX,nY,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: JAN 19, 2021
% LAST MODIFIED: JAN 19, 2021
% output coordinate matrix
% outcoordsMatrix(nX,nY,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  modX=1; modY=1;
  optargin = size(varargin,2);
  if(optargin==1)
    modX=varargin{1};
    modY=varargin{1};
  end
  if(optargin==2)
    modX=varargin{1};
    modY=varargin{2};
  end
  nX2=length(1:modX:nX); nY2=length(1:modY:nY);
  coords=zeros(nX2*nY2,2);
  cc=1;
  for i=1:modX:nX
    for j=1:modY:nY
      coords(cc,1)=i;
      coords(cc,2)=j;
      cc=cc+1;
    end
  end

