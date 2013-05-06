function f=cf(x,cflo,cfhi)
#  Cochlea to frequency functions for Cat model only
if nargin == 1
	cflo=200;
	cfhi=48000;
end

nchannels=100;
xlo = f2x(cflo);
xhi = f2x(cfhi);
xcenter = (xhi - xlo)/2.0 + xlo;
delx = (xhi-xlo)/(nchannels-1);
f = ceil(x2f(xlo + (x)*delx));
endfunction

function x=f2x(f)
x= 11.9 * log10(0.80 + (f / 456.0));
endfunction

function f = x2f(x)
f= 456.0*((10**(x/11.9))-0.80);
endfunction


