function make_anf_rate_vs(cflo,cfhi)

if nargin == 0
  cflo = 200;
  cfhi = 48000;
end

hsrfile = fopen("response_area-4.dat","w+");
lsrfile = fopen("response_area-5.dat","w+");
[status lsfiles] = system("find -maxdepth 1 -type d | grep -v -e '^.$' | tr -d './' |sort -n");
lsfiles = regexp(lsfiles,'\n','split');lsfiles(end)=[];

num_freqmod= numel(lsfiles);
SI=zeros(100,num_freqmod);  RT=zeros(100,num_freqmod);  
SI2=zeros(100,num_freqmod); RT2=zeros(100,num_freqmod); 
tic
for ii = 1:num_freqmod
  freqmod = str2num(lsfiles{1,ii})
  hsrras = load([lsfiles{1,ii} "/anHSR_raster.dat"] );
  lsrras = load([lsfiles{1,ii} "/anLSR_raster.dat"] );

      % rate= [];
      % for jj=0:99
      %   xh = hist(ras(find(ras(:,1)==jj),4),500);
      %  fprintf(rfile,"%d\t%d\t%d\t%g\n",freqmod,jj,cf(jj),mean(xh(100:end))*100/5);
      % end
      % fprintf(rfile,"\n");
      freqmod
      period = 1000/freqmod;
      for jj=0:99
	spk_index = find(hsrras(:,1)==jj);
	xsumSI = 0; ysumSI=0; n=0;
	for kk = 1:length(spk_index)
	  spike = hsrras(spk_index(kk),4);
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

	spk_index = find(lsrras(:,1)==jj);
	xsumSI = 0; ysumSI=0; n=0;
	for kk = 1:length(spk_index)
	  spike = lsrras(spk_index(kk),4);
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
            %% Rate 
      xh = hist(hsrras(find(hsrras(:,1)==jj),4),1:300);    
      fprintf(hsrfile,"%d\t%d\t%.3f\n",freqmod,cf(jj,cflo,cfhi),mean(xh(20:end))/0.28);
      %% LSR
      if numel(lsrras(find(lsrras(:,1)==jj)))>2
	xh = hist(lsrras(find(lsrras(:,1)==jj),4),1:300);
	%      LSR(ii,jj+1) = mean(xh(20:end-1))*100/5;
	fprintf(lsrfile,"%d\t%d\t%.3f\n",freqmod,cf(jj,cflo,cfhi),mean(xh(20:end))/0.28);
      else 
	fprintf(lsrfile,"%d\t%d\t%g\n",freqmod,cf(jj,cflo,cfhi),0.0);
      end
            
      end
      fprintf(hsrfile,"\n");
      fprintf(lsrfile,"\n");
end

fclose(hsrfile);
fclose(lsrfile);

    dlmwrite("vsSPIKES.4.dat",SI,"\t","precision","%.4f");  dlmwrite("rayltest.4.dat",RT,"\t","precision","%.2f");	  
    dlmwrite("vsSPIKES.5.dat",SI2,"\t","precision","%.4f"); dlmwrite("rayltest.5.dat",RT2,"\t","precision","%.2f");  
    clear SI SI2 RT RT2

toc

