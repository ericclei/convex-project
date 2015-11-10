function is = is_diagonally_dominant(A,strict)
%strict by default
if nargin < 2
  strict = 1;
end

if strict
  is = all( 2 * abs(diag(A)) > sum(abs(A),2) );
else
  is = all( 2 * abs(diag(A)) >= sum(abs(A),2) );
end