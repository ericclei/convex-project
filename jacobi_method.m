function [x,n_iter] = jacobi_method(A,b)
%[x,n_iter] = JACOBI_METHOD(A,b)
%solves Ax=b
%converges if A is diagonally dominant
n = size(A,1);
assert(size(A,2)==n);
D = spdiags(diag(A),0,n,n);
R = A-D;
x = D \ b;
tol = 1e-6;
n_iter = 0;
while 1
  n_iter = n_iter+1;
  xnew = D \ (b - R*x);
  if max(abs(xnew-x))<=tol
    x = xnew;
    break
  end
  x = xnew;
end