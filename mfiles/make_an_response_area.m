function make_an_response_area(cflo,cfhi)
%% make_an_response_area 2
% Make response area data for ANFs
%
% Michael Eager  (Mick.Eager@gmail.com)

if nargin == 0
  cflo = 200;
  cfhi = 48000;
end

onset=20;  % start measurement 20 ms after stim onset 

hsrfile = fopen("response_area.4.dat","w+");
lsrfile = fopen("response_area.5.dat","w+");
[status lsfiles] = system("find -maxdepth 1 -type d | grep -v -e '^.$' | tr -d './' |sort -n");
lsfiles = regexp(lsfiles,'\n','split');lsfiles(end)=[];
num_dir= numel(lsfiles);
tic

for ii = 1:num_dir
  dir_id = str2num(lsfiles{1,ii})
  hsrras = load([lsfiles{1,ii} "/anHSR_raster.dat"] );
  lsrras = load([lsfiles{1,ii} "/anLSR_raster.dat"] );
  num_hsr_fibres =   max(hsrras(:,2));
  num_lsr_fibres = max(lsrras(:,2));
  tstop = ceil(max(hsrras(:,4)));  % good chance the HSR units occur in last timestep, if not the max val

  dir_id
  for jj=0:99
    %% Rate 
    xh = hist(hsrras(find(hsrras(:,1)==jj),4),1:1:tstop);    
    fprintf(hsrfile,"%d\t%d\t%d\t%.3f\n",dir_id,jj,cf(jj,cflo,cfhi),mean(xh(20:end))/(num_hsr_fibres/1000));
    %% LSR
    if numel(lsrras(find(lsrras(:,1)==jj)))>2
      xh = hist(lsrras(find(lsrras(:,1)==jj),4),1:1:tstop);
      fprintf(lsrfile,"%d\t%d\t%d\t%.3f\n",dir_id,jj,cf(jj,cflo,cfhi),mean(xh(20:end))/(num_lsr_fibres/1000));
    else 
      fprintf(lsrfile,"%d\t%d\t%d\t%g\n",dir_id,jj,cf(jj,cflo,cfhi),0.0);
    end
    
    end
    fprintf(hsrfile,"\n");
    fprintf(lsrfile,"\n");
end

fclose(hsrfile);
fclose(lsrfile);

toc




