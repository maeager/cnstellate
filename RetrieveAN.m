% This m file can be used to extract AN filter bank responses to be used
% in plotting
ANfile = fopen(tmpANstr)
if (isopen(ANfile))
  tmp = fscanf(ANfile,"%d\n",6);       
  sg_rate = tmp(1);
  nsize = tmp(2);
  nchannels = tmp(3);
  cflo = tmp(4);
  cfhi = tmp(5);
  species = tmp(6);
  LSRout = zeroes(nchannels,nsize);
  HSRout = zeroes(nchannels,nsize);
  cf = zeroes(nchannels,1);
  for icf = 1:nchannels
    if (fscanf(ANfile,"%d\n",1) == icf)
      cf(icf) = fscanf(ANfile,"%d\n",1) 
      LSRout(:,icf) = fscanf(ANfile,"%.3f\t",nsize);
      HSRout(:,icf) = fscanf(ANfile,"%.3f\t",nsize);
    end
  end
  sg_tdres = 1/sg_rate
  fclose(ANfile)
end
