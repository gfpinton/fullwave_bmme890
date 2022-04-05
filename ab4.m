function [a b] = ab4(dx,kappax,alphax,ppw)
b=exp(-(dx./kappax+alphax).*6.0606e-09*12./ppw);
a=dx./(kappax.*(dx+kappax.*alphax)).*(b-1);

