function [Blam] = PlancksLaw(lambda,T)

h = 6.626e-34;
c = 299792458;
kB = 1.38065e-23;

Blam = (2.*h.*c^2./lambda.^5)./(exp(h.*c./(lambda.*kB.*T))-1);

end