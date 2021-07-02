function y = FixGaps(x)
% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.
%

% R. Pawlowicz 6/Nov/99
% MODIFIED by YK 7/May/2015

%interpolation
y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0; % ?


y(bd)=interp1(gd,x(gd),find(bd));

% Approximate first and last NaNs that can't be interpolated.
count=1;
if isnan(x(1))
   while isnan(x(count))
      count=count+1; 
   end
   y(1:count-1)=x(count);
end
count=length(x);
if isnan(x(end))
   while isnan(x(count))
      count=count-1; 
   end
   y(count+1:end)=x(count);
end

end
