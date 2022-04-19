function [pxducer1 pxducer2] = filtpxducer(pxducer,omega0,dT2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% FIRST WRITTEN: 2022-04-19
% LAST MODIFIED: 2022-04-19
% Filter harmonic component in pxducer
% Make sure transmit pulse has enough spectral separation/non overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nT2=size(pxducer,1);
  
  f0=omega0/2/pi;
  f=(0:nT2-1)/(nT2-1)/(dT2);

  flt=exp(-((f-f0)/f0*3).^4)';
  flt=exp(-((f-f0)/f0*3).^2)';
  flmat=flt*ones(1,size(pxducer,2));
  flt=exp(-((f-2*f0)/f0*3).^2)';
 %plot(f,flt), xlim([0 4*f0]), grid on, hold on
 % plot(f,dbzero(flt)), xlim([0 4*f0]), ylim([-60 0]),grid on, hold on
  flmat_h=flt*ones(1,size(pxducer,2));
    
  fpxducer=fft(double(pxducer));
  for k=1:size(pxducer,3)
    pxducer1(:,:,k)=single(ifft(fpxducer(:,:,k).*flmat,'symmetric'));
    pxducer2(:,:,k)=single(ifft(fpxducer(:,:,k).*flmat_h,'symmetric'));
  end
