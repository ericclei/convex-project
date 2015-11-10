function [x,n_iter] = successive_over_relaxation(A,b,w,verbose)
%[x,n_iter] = SUCCESSIVE_OVER_RELAXATION(A,b,w)
%solves Ax=b
%w is the relaxation parameter. w=1 is Gauss-Seidel
%if A==A' and 0<w<2, then SOR converges
if nargin < 4
  verbose = 0;
end
n = size(A,1);
assert(size(A,2)==n);
L = tril(A);
U = A - L;
x = L \ b;
tol = 1e-6;
n_iter = 0;
while 1
  n_iter = n_iter+1;
  x_gs = L \ (b - U*x);
  xnew = (1-w)*x+w*x_gs;
  if verbose && mod(n_iter,1000)==0
    disp(max(abs(xnew-x)));
  end
  if max(abs(xnew-x))<=tol
    x = xnew;
    break
  end
  x = xnew;
end