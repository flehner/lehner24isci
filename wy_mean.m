function water_year_mean=wy_mean(x)
% calculates water year weighted mean from monthly values
% assumptions: time is multiple of 12, starts at January, year has 365 days
% --> not implemented yet: multidimensional input matrix (e.g., a time, lat, lon field)

wy_w = [31 30 31 31 28 31 30 31 30 31 31 30]/365;

xs = size(x);
if xs(1) > xs(2)
  x = x';
end

ntime = floor(length(x)/12)-1;
water_year_mean = NaN(1,ntime);
for i = 1:ntime
  water_year_mean(i+1) = sum(x((i-1)*12+10:i*12+9).*wy_w);
end
water_year_mean = water_year_mean(2:end);

return
