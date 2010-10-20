function getnotchcheb(Fs,fstop1,notchwidth,Rs,filename,printflag)
%close all;
%Fs = 100e3
x = 2e-5*10^(100/20)*randn(1,Fs/10); %get 100 ms of noise from Fs=100KHz
xstep=ceil(2*Fs/1000);    % one spectral slice every 20 ms
window=ceil(10*Fs/1000); % 10 ms data window
if nargin == 6 && printflag == 1
  subplot(3,1,1);
  specgram(x, 2^nextpow2(window), Fs, window, window-xstep);
end
 
%notchwidth=1; n=9
%notchwidth=1/2; n=9
%notchwidth=1/4; n=7
%Rs= 30; % dB
if notchwidth > 10 
  printf("Two stopbands mode\n")
  if notchwidth > fstop1
    fstop2=notchwidth; 
  else
    printf("Error Two stopbands mode, fstop2 not bigger than fstop1.\n")
  end 
else
  fstop2=fstop1+notchwidth*fstop1; 
end
% fstop1 = fstop2-notchwidth*fstop2		
Ws1 = fstop1/(Fs/2);
Ws2 = fstop2/(Fs/2);
if Ws1 < 0
   Ws1=0;
end
if Ws2 > 1
  Ws2=1;
end

Wp1 = (100+fstop1)/(Fs/2);
Wp2 = (fstop2-10)/(Fs/2);
Rp=0.5;
% Find the best order for filter
%[n, Wc] = cheb2ord([Wp1, Wp2], [Ws1, Ws2], Rp, Rs)
% Band pass (Ws1<Wp1<Wp2<Ws2) or band reject (Wp1<Ws1<Ws2<Wp2)
%filter design. Wp gives the edges of the pass band, and Ws gives
%the edges of the stop band.
% [n, Wc] = cheb2ord (1000/(Fs/2), 1200/(Fs/2), 0.5, Rs) ;
n=14
do
printf("Ws1 %g  Ws2 %g\n",Ws1,Ws2)
[b, a] = cheby2(n, Rs, [Ws1, Ws2], 'stop');
% band reject filter with edges pi*Wl and pi*Wh radians
%[b, a] = cheby2 (n, 29, Wc);
[h, w] = freqz (b, a, [], Fs);
n-=1
until (max(20*log10(abs(h))) < 0.5)
notchnoise = filter(b,a,x);

if nargin == 6 && printflag == 1  
  subplot(3,1,2);
  plot (w, 20*log10(abs(h)), ";;");
  ylabel("Attenuation (dB)","fontsize",14);
  xlabel("Frequency (Hz)","fontsize",14);
  subplot(3,1,3)
  specgram(notchnoise, 2^nextpow2(window), Fs, window, window-xstep);
end
%save(filename, "notchnoise");
f = fopen(filename,"w+")
fprintf(f,"%d\n%d\n",1,length(notchnoise))
for ii=1:length(notchnoise)
    fprintf(f,"%.4f\t",notchnoise(ii))
end
fclose(f)
if nargin == 6 && printflag == 1  
  print("-depsc2",[filename ".eps"]);
end


