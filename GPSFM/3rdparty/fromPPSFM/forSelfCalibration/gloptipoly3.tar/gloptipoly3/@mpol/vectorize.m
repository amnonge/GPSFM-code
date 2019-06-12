function str = vectorize(p,varx,vary,varz)
% @MPOL/VECTORIZE Vectorize polynomial
%
% IF P is a univariate, bivariate or trivariate polynomial (class MPOL), the instruction
%
%   VECTORIZE(P) 
%
% returns a character string (class CHAR) corresponding to the polynomial
% with variables X (and Y and Z). If X (and Y and Z) are existing Matlab vectors
% or matrices (class DOUBLE) then EVAL(VECTORIZE(P)) returns the evaluation
% of P at the values of X (and Y and Z).
%
% With the syntax VECTORIZE(P,VARX,VARY,VARZ) the names of the
% variables can be specified in strings VARX, VARY and VARZ.

% D. Henrion, 6 April 2010
% Last modified on 24 January 2014

if nargin < 2
 varx = 'x';
end
if nargin < 3
 vary = 'y';
end
if nargin < 4
 varz = 'z';
end

N = length(p);
if N > 1
 error('Scalar polynomials only')
end
n = length(indvar(p));
if n > 3
 error('Univariate, bivariate or trivariate polynomials only')
end
pp = pow(p); cp = coef(p);
[q,m] = size(pp);
str = '';
for k = 1:q
 str = [str '+(' num2str(cp(k)) ').*' varx '.^(' int2str(pp(k,1)) ')'];
 if m > 1
  str = [str '.*' vary '.^(' int2str(pp(k,2)) ')'];
 end
 if m > 2
  str = [str '.*' varz '.^(' int2str(pp(k,3)) ')'];
 end
end 


 




