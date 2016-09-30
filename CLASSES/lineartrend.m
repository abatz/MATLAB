function [TREND,tval]=lineartrend(DATA,pvalue);
% function [TREND,tval]=lineartrend(DATA,pvalue);
% calculates TREND and significance of a time series.  Trend is reported
% in units per time unit, significance = 1 if trend is significant, 0 if
% not significant at p<pvalue 
     d1=squeeze(DATA);
     D4=1:length(d1);D4=D4';
     rtt=0.75;  % only calculate trend is at least 75% of data are valid
     if length(a11)>=rtt*length(DATA)
     c=polyfit(D4(a11),d1(a11),1);
     TREND=c(1);  % trend per time unit
     % estimate autocorrelation of residuals and statistical significance
     refit=polyval(c,D4(a11));
     errors=refit-d1(a11);
     se=sqrt(1/(length(a11)-2)*sum(errors(:).^2));
     meant=mean(D4(a11));
     terrors=(D4(a11))-meant;
     sb=se/sqrt(sum(terrors(:).^2));
     autor=auto(errors,1);autor=autor(2);
     neff=ceil((length(a11)-1+1)*(1-autor)/(1+autor));
     tval=tinv(1-pvalue/2,neff)*sb;
     else TREND=NaN;tval=0;end
if abs(TREND)>tval tval=1;else tval=0;end


% This function calculates the autocorrelation function for the 
% data set "origdata" for a lag timescale of 0 to "endlag" and outputs 
% the autocorrelation function in to "a".
%function a = auto( origdata, endlag);

function a = auto(data, endlag);

N=length(data);

%now solve for autocorrelation for time lags from zero to endlag
for lag=0:endlag
  data1=data(1:N-lag);
   data1=data1-mean(data1);
  data2=data(1+lag:N);
   data2=data2-mean(data2);
  a(lag+1) = sum(data1.*data2)/sqrt(sum(data1.^2).*sum(data2.^2));
  clear data1
  clear data2
end
  

function x = tinv(p,v);
%TINV   Inverse of Student's T cumulative distribution function (cdf).
%   X=TINV(P,V) returns the inverse of Student's T cdf with V degrees 
%   of freedom, at the values in P.
%
%   The size of X is the common size of P and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    


if nargin < 2, 
    error('Requires two input arguments.'); 
end

[errorcode p v] = distchck(2,p,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
x=zeros(size(p));

k = find(v < 0  | v ~= round(v));
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end

k = find(v == 1);
if any(k)
  x(k) = tan(pi * (p(k) - 0.5));
end

% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
k0 = find(p == 0);
if any(k0)
    tmp   = Inf;
    x(k0) = -tmp(ones(size(k0)));
end
k1 = find(p ==1);
if any(k1)
    tmp   = Inf;
    x(k1) = tmp(ones(size(k1)));
end

k = find(p >= 0.5 & p < 1);
if any(k)
    z = betainv(2*(1-p(k)),v(k)/2,0.5);
    x(k) = sqrt(v(k) ./ z - v(k));
end

k = find(p < 0.5 & p > 0);
if any(k)
    z = betainv(2*(p(k)),v(k)/2,0.5);
    x(k) = -sqrt(v(k) ./ z - v(k));
end
function x = betainv(p,a,b);
%BETAINV Inverse of the beta cumulative distribution function (cdf).
%	X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%	parameters A and B at the values in P.
%
%	The size of X is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 
%
%	BETAINV uses Newton's method to converge to the solution.

%	Reference:
%	   [1]     M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964.

%	B.A. Jones 1-12-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:53:18 $

if nargin <  3, 
    error('Requires three input arguments.'); 
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

%   Initialize x to zero.
x = zeros(size(p));

%   Return NaN if the arguments are outside their respective limits.
k = find(p < 0 | p > 1 | a <= 0 | b <= 0);
if any(k),
   x(k) = NaN * ones(size(k)); 
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0), 
    x(k0) = zeros(size(k0)); 
end

k1 = find(p==1);
if any(k1), 
    x(k1) = ones(size(k1)); 
end

% Newton's Method.
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;

k = find(p > 0 & p < 1 & a > 0 & b > 0);
pk = p(k);

%   Use the mean as a starting guess. 
xk = a(k) ./ (a(k) + b(k));


% Move starting values away from the boundaries.
if xk == 0,
    xk = sqrt(eps);
end
if xk == 1,
    xk = 1 - sqrt(eps);
end


h = ones(size(pk));
crit = sqrt(eps); 

% Break out of the iteration loop for the following:
%  1) The last update is very small (compared to x).
%  2) The last update is very small (compared to 100*eps).
%  3) There are more than 100 iterations. This should NEVER happen. 

while(any(abs(h) > crit * abs(xk)) & max(abs(h)) > crit    ...
                                 & count < count_limit), 
                                 
    count = count+1;    
    h = (betacdf(xk,a(k),b(k)) - pk) ./ betapdf(xk,a(k),b(k));
    xnew = xk - h;

% Make sure that the values stay inside the bounds.
% Initially, Newton's Method may take big steps.
    ksmall = find(xnew < 0);
    klarge = find(xnew > 1);
    if any(ksmall) | any(klarge)
        xnew(ksmall) = xk(ksmall) /10;
        xnew(klarge) = 1 - (1 - xk(klarge))/10;
    end

    xk = xnew;  
end

% Return the converged value(s).
x(k) = xk;

if count==count_limit, 
    fprintf('\nWarning: BETAINV did not converge.\n');
    fprintf('The last Newton step was:\n');
    disp(h);
end

function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
%DISTCHCK Checks the argument list for the probability functions.

%	B.A. Jones  1-22-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:54:03 $

errorcode = 0;

if nparms == 1
    out1 = arg1;
    return;
end
    
if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end     
    end
    if scalararg1
        out1 = arg1(ones(r2,1),ones(c2,1));
    else
        out1 = arg1;
    end
    if scalararg2
        out2 = arg2(ones(r1,1),ones(c1,1));
    else
        out2 = arg2;
    end
end
    
if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1
        [rows columns] = size(arg1);
    elseif ~scalararg2
    [rows columns] = size(arg2);
    else
        [rows columns] = size(arg3);
    end
    out1 = arg1(ones(rows,1),ones(columns,1));
    out2 = arg2(ones(rows,1),ones(columns,1));
    out3 = arg3(ones(rows,1),ones(columns,1));
    out4 =[];
end

if nparms == 4
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    [r4 c4] = size(arg4);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);
    scalararg4 = (prod(size(arg4)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg1 & ~scalararg4
        if r1 ~= r4 | c1 ~= c4
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg4 & ~scalararg2
        if r4 ~= r2 | c4 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg3 & ~scalararg4
        if r3 ~= r4 | c3 ~= c4
            errorcode = 1;
            return;         
        end
    end


    if ~scalararg1
        [rows columns] = size(arg1);
    elseif ~scalararg2
    [rows columns] = size(arg2);
    elseif ~scalararg3
        [rows columns] = size(arg3);
    else
        [rows columns] = size(arg4);
    end
    out1 = arg1(ones(rows,1),ones(columns,1));
    out2 = arg2(ones(rows,1),ones(columns,1));
    out3 = arg3(ones(rows,1),ones(columns,1));
    out4 = arg4(ones(rows,1),ones(columns,1));
end

function p = betacdf(x,a,b);
%BETACDF Beta cumulative distribution function.
%	P = BETACDF(X,A,B) returns the beta cumulative distribution
%	function with parameters A and B at the values in X.
%
%	The size of P is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 
%
%	BETAINC does the computational work.

%	Reference:
%	   [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964, 26.5.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:53:16 $

if nargin<3, 
   error('Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
   error('The arguments must be the same size or be scalars.');
end

% Initialize P to 0.
p = zeros(size(x));

k1 = find(a<=0 | b<=0);
if any(k1)
   p(k1) = NaN * ones(size(k1)); 
end

% If is X >= 1 the cdf of X is 1. 
k2 = find(x >= 1);
if any(k2)
   p(k2) = ones(size(k2));
end

k = find(x > 0 & x < 1 & a > 0 & b > 0);
if any(k)
   p(k) = betainc(x(k),a(k),b(k));
end

% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
p(k) = ones(size(k));

function y = betapdf(x,a,b)
%BETAPDF Beta probability density function.
%	Y = BETAPDF(X,A,B) returns the beta probability density 
%	function with parameters A and B at the values in X.
%
%	The size of Y is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 

%	References:
%	   [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964, 26.1.33.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:53:21 $

if nargin < 3, 
   error('Requires three input arguments.');
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y = zeros(size(x));

% Return NaN for parameter values outside their respective limits.
k1 = find(a <= 0 | b <= 0 | x < 0 | x > 1);
if any(k1)
    y(k1) = NaN * ones(size(k1)); 
end

% Return the beta density function for valid parameters.
k = find(~(a <= 0 | b <= 0 | x < 0 | x > 1));
if any(k)
    y(k) = x(k) .^ (a(k) - 1) .* (1 - x(k)) .^ (b(k) - 1) ./ beta(a(k),b(k));
end

