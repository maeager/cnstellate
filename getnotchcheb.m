function getnotchcheb(Fs,fstop1,notchwidth,Rs,filename)
%close all;
%Fs = 100e3
x = 2e-5*10^(100/20)*randn(1,Fs/10);
xstep=ceil(2*Fs/1000);    # one spectral slice every 20 ms
window=ceil(10*Fs/1000); # 10 ms data window
subplot(3,1,1);
specgram(x, 2^nextpow2(window), Fs, window, window-xstep);
n=9

%notchwidth=1; n=9
%notchwidth=1/2; n=9
%notchwidth=1/4; n=7
%Rs= 30; % dB
fstop2=fstop1+notchwidth*fstop1; 
%% fstop1 = fstop2-notchwidth*fstop2		
		
Ws1 = fstop1/(Fs/2)
Ws2 = fstop2/(Fs/2)
Wp1 = 0
Wp2 = 1
Rp=0.5

% Find the best order for filter
%[n, Wc] = cheb2ord([Wp1, Wp2], [Ws1, Ws2], Rp, Rs)
% Band pass (Ws1<Wp1<Wp2<Ws2) or band reject (Wp1<Ws1<Ws2<Wp2)
%filter design. Wp gives the edges of the pass band, and Ws gives
%the edges of the stop band.
% [n, Wc] = cheb2ord (1000/(Fs/2), 1200/(Fs/2), 0.5, Rs);
[b, a] = cheby2(n, Rs, [Ws1, Ws2], 'stop');
  # band reject filter with edges pi*Wl and pi*Wh radians
%[b, a] = cheby2 (n, 29, Wc);
subplot(3,1,2)
 [h, w] = freqz (b, a, [], Fs);
 plot (w, 20*log10(abs(h)), ";;");
ylabel("Attenuation (dB)","fontsize",14);
 xlabel("Frequency (Hz)","fontsize",14);

subplot(3,1,3)
notchnoise = filter(b,a,x);
specgram(notchnoise, 2^nextpow2(window), Fs, window, window-xstep);

%save(filename, "notchnoise");
f = fopen(filename,"w+")
for ii=1:length(notchnoise)
    fprintf(f,"%.4f\t",notchnoise(ii))
end
fclose(f)
print("-depsc2",[filename ".eps"]);
