function y = rm(x,wl,dim);
% Computes honest running mean without any extrapolation - where mean can't be
% computed due to window length there is NaN
% First value is plotted at rounded-down halfway point of selected window length
% - x   = data (1- or 2-dimensional)
% - wl  = window length for each mean
% - dim = array dimension over which to compute running mean
%         (optinal argument; only needed if 2-dimensional)
if wl == 1
  % -- do nothing
  if nargin == 2
    y = x';
  else
    y = x;
  end
else
  if nargin == 2
    y     = NaN(1,length(x));
    for i = 1:length(x)-wl+1
      y(i+round(wl/2)) = nanmean(x(i:i+wl-1));
    end
  elseif nargin == 3
    ndim  = size(x);
    y     = NaN(ndim);
    if dim == 1
      for j = 1:ndim(2)
        for i = 1:ndim(1)-wl+1
          y(i+round(wl/2),j) = nanmean(x(i:i+wl-1,j),dim);
        end
      end
    else % if dim == 2
      for j = 1:ndim(2)-wl+1
        for i = 1:ndim(1)
          y(i,j+round(wl/2)) = nanmean(x(i,j:j+wl-1),dim);
        end
      end
    end
  end
end


return
