
rfile = fopen("response_area.4.dat","w+")
HSR_SI = zeros(100,16);
for ii = 1:16
  freqmod = (ii-1)*50 + 50
  ras = load([num2str(freqmod) "/anHSR_raster.dat"] );

 % rate= [];
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
clear HSR_SI

LSR_SI = zeros(100,16);
rfile = fopen("response_area.5.dat","w+")
for ii = 1:16
  freqmod = (ii-1)*50 + 50
  ras = load([num2str(freqmod) "/anLSR_raster.dat"] );

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
   LSR_SI(jj+1,ii) = sqrt(xsumSI^2 + ysumSI^2)/n;
 end
end


fclose(rfile);
dlmwrite("vsSPIKES.5.dat",LSR_SI,"\t");




