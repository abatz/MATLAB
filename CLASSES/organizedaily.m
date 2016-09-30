function [dataout]=organizedata(day,month,year,data);
% function [dataout]=organizedata(day,month,year,data);
% organizes data into matrix of data organized by day of year and year
doy=dayofyear_fixed(year,month,day);
dataout=NaN*ones(366,max(year)-min(year)+1);
minyear=min(year);
for i=1:366
 f=find(doy==i);
 dataout(i,year(f)+1-minyear)=data(f);
end


function yd = dayofyear_fixed(varargin)
%DAYOFYEAR Ordinal number of day in a year.
%
%   DAYOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
%   day number in the given year plus a fractional part depending on the
%   time of day.
%
%   Any missing MONTH or DAY will be replaced by 1.  HOUR, MINUTE or SECOND
%   will be replaced by zeros.
%
%   If no date is specified, the current date and time is used.  Gregorian
%   calendar is assumed.


   nargsin = nargin;
   error(nargchk(0, 6, nargsin));
   if nargsin
      argv = { 1 1 1 0 0 0 };
      argv(1:nargsin) = varargin;
   else
      argv = num2cell(clock);
   end
   [year, month, day, hour, minute, second] = deal(argv{:});
   year=year(:);
   month=month(:);
   day=day(:);

   days_in_prev_months = [0 31 60 91 121 152 182 213 244 274 305 335];
   yd = days_in_prev_months(month)' ...               % days in prev. months
        + day ...                                    % day in month
        + ( second + 60*minute + 3600*hour )/86400;  % part of day
   yd=yd;%+365;%*(year-1997);





