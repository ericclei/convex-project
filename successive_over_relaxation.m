function [x,n_iter,err_path] = successive_over_relaxation(A,b,w,verbose,x_star)
%[x,n_iter] = SUCCESSIVE_OVER_RELAXATION(A,b,w)
%solves Ax=b
%w is the relaxation parameter. w=1 is Gauss-Seidel
%if A==A' and 0<w<2, then SOR converges
if ~exist('verbose','var')
  verbose = 0;
end
have_x_star = exist('x_star','var');
n = size(A,1);
assert(size(A,2)==n);
if have_x_star, err_path = nan(1e9,1); end
L = tril(A);
U = A - L;
invL = L^-1;
x = invL * b;
if have_x_star, err_path(1) = norm(x-x_star); end
tol = 1e-6;
n_iter = 0;
while 1
  n_iter = n_iter+1;
  z = b-U*x;
  x_gs = invL*z;
  xnew = (1-w)*x+w*x_gs;
  if have_x_star, err_path(1+n_iter) = norm(xnew-x_star); end
  if verbose && mod(n_iter,10000)==0
    disp(n_iter);
    disp(max(abs(xnew-x)));
  end
  if max(abs(xnew-x))<=tol
    x = xnew;
    break
  end
  x = xnew;
end
if have_x_star, err_path = err_path(1:n_iter+1); end