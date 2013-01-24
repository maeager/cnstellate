function getbandcheb(Fs,fpass1,bandwidth,Rs,filename,printflag,n_)
%close all;
%Fs = 100e3
x = 2e-5*10^(100/20)*randn(1,Fs/10); %get 100 ms of noise from Fs=100KHz
xstep=ceil(2*Fs/1000);    % one spectral slice every 20 ms
window=ceil(10*Fs/1000); % 10 ms data window
if nargin == 6 && printflag == 1
  subplot(3,1,1);
  specgram(x, 2^nextpow2(window), Fs, window, window-xstep);
end
 
%bandwidth=1; n=9
%bandwidth=1/2; n=9
%bandwidth=1/4; n=7
%Rs= 30; % dB
if bandwidth > 10 
  printf("Two passbands mode\n")
  if bandwidth > fpass1
    fpass2=bandwidth; 
  else
    printf("Error Two passbands mode, fpass2 not bigger than fpass1.\n")
  end 
else
  fpass2=fpass1+bandwidth*fpass1; 
end
% fpass1 = fpass2-bandwidth*fpass2		
Wp1 = fpass1/(Fs/2);
Wp2 = fpass2/(Fs/2);
if Wp1 < 0
   Wp1=0;
end
if Wp2 > 1
  Wp2=1;
end

Ws1=Wp1-0.01;
Ws2=Wp2+0.01;
%Ws1=0;%
%Ws2=1;%

Rp=0.5;
% Find the best order for filter
if nargin ~= 7
   [n, Wc] = cheb2ord([Wp1, Wp2], [Ws1, Ws2], Rp, Rs)
else
    n = n_;
end 
% Band pass (Ws1<Wp1<Wp2<Ws2) or band reject (Wp1<Ws1<Ws2<Wp2)
%filter design. Wp gives the edges of the pass band, and Ws gives
%the edges of the pass band.
% [n, Wc] = cheb2ord (1000/(Fs/2), 1200/(Fs/2), 0.5, Rs) ;
printf("Ws1 %g  Ws2 %g\n",Ws1,Ws2)
[b, a] = cheby2(n, Rs, [Ws1, Ws2]);

[h, w] = freqz (b, a, [], Fs);
bandnoise = filter(b,a,x);

if nargin == 6 && printflag == 1  
  subplot(3,1,2);
  plot (w, 20*log10(abs(h)), ";;");
  ylabel("Attenuation (dB)","fontsize",14);
  xlabel("Frequency (Hz)","fontsize",14);
  subplot(3,1,3)
  specgram(bandnoise, 2^nextpow2(window), Fs, window, window-xstep);
end
%save(filename, "bandnoise");
f = fopen(filename,"w+")
fprintf(f,"%d\n%d\n",1,length(bandnoise))
for ii=1:length(bandnoise)
    fprintf(f,"%.4f\t",bandnoise(ii))
end
fclose(f)
if nargin == 6 && printflag == 1  
  print("-depsc2",[filename ".eps"]);
end


