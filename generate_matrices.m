clear
home
close all
%% generate a sparse ill-conditioned positive definite matrix
n = 1e5;
density = 2/n;
cond = 1e4/eps;
rc = 1/cond;
m = 10;
Matrices = cell(m,1);
%%
for i=2:m
  tic
  Matrices{i} = sprandsym(n, density, rc, 1);
  toc
end
% nnzwanted = round(density*n*n);
% ratio = rc ^ (1/(n-1));
% R = sparse(1:n,1:n, (ratio .^ ( 0:(n-1)))  );
% nnzr = n;
% %   random jacobi rotations
% while ( nnzr < .95*nnzwanted )
%     R = rjr(R);
%     nnzr = nnz(R);
% end
%%
save matrices