


## function [ampl]=dB2ampl(dB)
## 	ampl = 10^(dB/20)
## endfunction
	
## function [ampl]=dBSPL2ampl(dB)
## 	ampl = 2e-5*10^(dB/20)
## endfunction


%function getnotch()
close all

Fs = 1e5;
notchwidth=1/4;
notchdepth= -30; % dB
depth= dB2ampl(notchdepth);
fstop1=12.5e3;
fstop2 = fstop1+notchwidth*fstop1
sb1 = fstop1/(Fs/2)
sb2 = fstop2/(Fs/2)

f = [0 sb1 sb1 sb2 sb2   1]
m = [1 1 depth depth 1 1]
[b,a] = fir2(100,f,m)
  [h, w] = freqz(b,a);
  plot(f,m,';target response;',w/pi,abs(h),';filter response;');

