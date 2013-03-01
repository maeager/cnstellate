
% This m file can be used to extract AN filter bank responses to be used
% in plotting

% 
%   
%    Copyright © 2012 Michael Eager, (mick.eager@gmail.com)
%    cnstellate was written as part of my PhD at The University of Melbourne
%
%    This file is part of cnstellate.
% 
%    This is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.  
%


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
