function f0synchindex()


freqmod = 150
SI=zeros(100,15);  RT=10e9*ones(100,15); 
SI2=zeros(100,15); RT2=10e9*ones(100,15);

for spl = 0:5:70
  ras = load([num2str(spl) "/anHSR_raster.dat"] );
  ras2 = load([num2str(spl) "/anLSR_raster.dat"] );
 # rate= [];
# for jj=0:99
#   xh = hist(ras(find(ras(:,1)==jj),4),500);
#  fprintf(rfile,"%d\t%d\t%d\t%g\n",freqmod,jj,cf(jj),mean(xh(100:end))*100/5);
# end
# fprintf(rfile,"\n");

 for jj=0:99
   spk_index = find(ras(:,1)==jj);
   xsumSI = 0; ysumSI=0; n=0;
   for kk = 1:length(spk_index)
     spike = ras(spk_index(kk),4);
     if spike > 40
       theta = 2*pi* mod(spike , 1000/freqmod) / (1000/freqmod);
       xsumSI += cos(theta);
       ysumSI += sin(theta);
       n+=1;
     end
   end
   if n > 0
     SI(jj+1,spl/5 + 1) = sqrt(xsumSI^2 + ysumSI^2)/n;
     RT(jj+1,spl/5 + 1) = 2*n*SI(jj+1,spl/5 + 1)^2;
   end

   spk_index = find(ras2(:,1)==jj);
   xsumSI = 0; ysumSI=0; n=0;
   for kk = 1:length(spk_index)
     spike = ras2(spk_index(kk),4);
     if spike > 40
       theta = 2*pi* mod(spike , 1000/freqmod) / (1000/freqmod);
       xsumSI += cos(theta);
       ysumSI += sin(theta);
       n+=1;
     end
   end
   if n > 0
     SI2(jj+1,spl/5 + 1) = sqrt(xsumSI^2 + ysumSI^2)/n;
     RT2(jj+1,spl/5 + 1) = 2*n*SI2(jj+1,spl/5 + 1)^2;
   end
 end
end


dlmwrite("vsSPIKES.4.dat",SI,"\t");  dlmwrite("rayltest.4.dat",RT,"\t");  
dlmwrite("vsSPIKES.5.dat",SI2,"\t"); dlmwrite("rayltest.5.dat",RT2,"\t"); 
clear SI SI2 RT RT2

binwidth = 0.1  % ms
period = 1000 / freqmod %ms

for ii=0:3
  SI=zeros(100,15); RT=10e9*ones(100,15);
  for spl = 0:5:70
    phist = load([num2str(spl) "/periodhist." num2str(ii) ".dat"] );

# rate= [];
# for jj=0:99
#   xh = hist(ras(find(ras(:,1)==jj),4),500);
#  fprintf(rfile,"%d\t%d\t%d\t%g\n",freqmod,jj,cf(jj),mean(xh(100:end))*100/5);
# end
# fprintf(rfile,"\n");

    for jj=0:99
      spk_index = find(phist(:,1)==jj);
      xsumSI = 0; ysumSI=0; n=0;
      for kk = 1:length(spk_index)
	phasetime = phist(spk_index(kk),2);
 
	theta = 2*pi* mod(phasetime, period) / period;
	xsumSI += cos(theta)* phist(spk_index(kk),3);
	ysumSI += sin(theta)* phist(spk_index(kk),3);
      end
      n= sum(phist(spk_index(:),3));
      if n > 0 
	SI(jj+1,spl/5 + 1) = sqrt(xsumSI^2 + ysumSI^2)/n;
	RT(jj+1,spl/5 + 1) = 2*n*SI(jj+1,spl/5 + 1)^2;
      end

    end
  end

dlmwrite(["vsSPIKES." num2str(ii) ".dat"],SI,"\t");
dlmwrite(["rayltest." num2str(ii) ".dat"],RT,"\t");

end

endfunction

