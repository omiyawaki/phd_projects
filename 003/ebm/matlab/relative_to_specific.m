function q = relative_to_specific(h,T,p)

% This function converts relative humidities to specific humidities via a
% formula derived in my 3/23/17 notes.  Input variables h, T, and p can
% have any dimensions but must be able to be multiplied and added as needed
% via bsxfun.  In other words, each dimension must be the same length in
% all three variables or a dimension can have a nonzero length in exactly
% one variable as long as it has length zero for the other two variables.
%

% INPUTS:
% h     Relative humidity, units [] (values should range from 0 to 1)
% T     Temperature, units K
% p     Total air pressure, units Pa
%
% OUTPUT:
% q     Specific humidity, units []

% 1) LOAD CONSTANTS
R         = todd_constants('R');
Rv        = todd_constants('Rv');
Rv_over_R = Rv/R;

% 2) GET SATURATION VAPOR PRESSURES
es        = compute_es(T);

% 3) COMPUTE SPECIFIC HUMIDITY
q = 1./(1+(Rv_over_R*(bsxfun(@rdivide,p,bsxfun(@times,h,es))-1)));