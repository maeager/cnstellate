
rfile = fopen("response_area.4.dat","w+")
for ii = 1:16
  freqmod = (ii-1)*50 + 50
  ras = load([num2str(freqmod) "/anHSR_raster.dat"] );

 # rate= [];
 for jj=0:99
   xh = hist(ras(find(ras(:,1)==jj),4),500);
  fprintf(rfile,"%d\t%d\t%d\t%g\n",freqmod,jj,cf(jj),mean(xh(100:end))*100/5);
 end
 fprintf(rfile,"\n");

 for jj=0:99
   spk_index = find(ras(:,1)==jj);
   xsumSI = 0; ysumSI=0; n=0;
   for kk = 1:length(spk_index)
     spike = ras(spk_index(kk),4);
     if spike > 40
       theta = 2*pi* mod(spike , 1000/freqmod);
       xsumSI += cos(theta);
       ysumSI += sin(theta);
       n+=1;
     end
   end
   HSR_SI(jj+1,ii) = sqrt(xsumSI^2 + ysumSI^2)/n;
 end
end


fclose(rfile);
dlmwrite("vsSPIKES.4.dat",HSR_SI,"\t");


#  Cochlea to frequency functions for Cat model only

function f=cf(x)
# Cat data
cflo=200
cfhi=48000
nchannels=100
xlo = f2x(cflo)
xhi = f2x(cfhi)
xcenter = (xhi - xlo)/2.0 + xlo
delx = (xhi-xlo)/(nchannels-1)
f = ceil(x2f(xlo + (x)*delx))
endfunction

function x=f2x(f)  
x= 11.9 * log10(0.80 + (f / 456.0))
endfunction

function f = x2f(x) 
f= 456.0*((10**(x/11.9))-0.80)
endfunction
