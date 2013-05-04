function make_an_vsSPIKES()

  for spl = 20:20:80

    SI=zeros(100,17);  RT=zeros(100,17);  
    SI2=zeros(100,17); RT2=zeros(100,17); 
    for freqmod = 50:50:1200

      ras = load([num2str(spl) "/" num2str(freqmod) "/anHSR_raster.dat"] );
      ras2 = load([num2str(spl) "/" num2str(freqmod) "/anLSR_raster.dat"] );
      # rate= [];
      # for jj=0:99
      #   xh = hist(ras(find(ras(:,1)==jj),4),500);
      #  fprintf(rfile,"%d\t%d\t%d\t%g\n",freqmod,jj,cf(jj),mean(xh(100:end))*100/5);
      # end
      # fprintf(rfile,"\n");
      freqmod;
      period = 1000/freqmod;
      for jj=0:99
	spk_index = find(ras(:,1)==jj);
	xsumSI = 0; ysumSI=0; n=0;
	for kk = 1:length(spk_index)
	  spike = ras(spk_index(kk),4);
	  if spike > 40
	    theta = 2*pi* mod(spike , period) / (period);
	    xsumSI += cos(theta);
	    ysumSI += sin(theta);
            n+=1;
	  end
	end
	if n > 4
	  SI(jj+1,freqmod/50) = sqrt(xsumSI^2 + ysumSI^2)/n;
	  RT(jj+1,freqmod/50) =  2*n*SI(jj+1,freqmod/50)^2;
	end

	spk_index = find(ras2(:,1)==jj);
	xsumSI = 0; ysumSI=0; n=0;
	for kk = 1:length(spk_index)
	  spike = ras2(spk_index(kk),4);
	  if spike > 40
	    theta = 2*pi* mod(spike , period) / (period);
	    xsumSI += cos(theta);
	    ysumSI += sin(theta);
	    n+=1;
	  end
	end
	if n > 4
	  SI2(jj+1,freqmod/50) = sqrt(xsumSI^2 + ysumSI^2)/n;
	  RT2(jj+1,freqmod/50) =  2*n*SI(jj+1,freqmod/50)^2;
	end
      end
    end


    dlmwrite([num2str(spl) "/vsSPIKES.4.dat"],SI,"\t");  dlmwrite([num2str(spl) "/rayltest.4.dat"],RT,"\t");	  
    dlmwrite([num2str(spl) "/vsSPIKES.5.dat"],SI2,"\t"); dlmwrite([num2str(spl) "/rayltest.5.dat"],RT2,"\t");  
    clear SI SI2 RT RT2

  end
