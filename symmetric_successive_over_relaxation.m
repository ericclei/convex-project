function [x,n_iter,err_path] = symmetric_successive_over_relaxation(A,b,w,verbose,x_star)

if ~exist('verbose','var')
  verbose = 0;
end
have_x_star = exist('x_star','var');
if nargout==3 && ~have_x_star
  error('need x_star to give err_path')
end
tol = 1e-5;
n = size(A,1);
assert(size(A,2)==n);
if have_x_star, err_path = nan(1e9,1); end
D = spdiags(diag(A),0,n,n);
L = tril(A-D);
L = sparse(L);
U = A-D-L;
U = sparse(U);
invDwL = (D+w*L)^-1;
x = invDwL*w*b;
if have_x_star
  norm_x_star = norm(x_star);
  err_path(1) = norm(x-x_star)/norm_x_star; 
end
n_iter = 0;
while 1
  n_iter = n_iter+1;
  if w==1
    xnew = invDwL*(b-U*x);
  else
    xnew = invDwL*(w*b-(w*U+(w-1)*D)*x);
  end
  if have_x_star, err_path(1+n_iter) = norm(xnew-x_star)/norm_x_star; end
  delta = norm(xnew-x)/norm(x);
  del(n_iter) = delta;
  if verbose && mod(n_iter,1000)==0
    disp(n_iter);
    disp(delta);
  end
  if delta<=tol
    x = xnew;
    break
  end
  x = xnew;
end
if have_x_star, err_path = err_path(1:n_iter+1); end