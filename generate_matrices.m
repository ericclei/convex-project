clear
home
close all
%% generate a sparse ill-conditioned positive definite matrix
n = 1e4;
density = 2/n;
cond = 2/eps;
rc = 1/cond;
m = 10;
Matrices = cell(m,1);
Bs = cell(m,1);
%%
for i=1:m
  tic
  A = sprandsym(n, density, rc, 1);
  Matrices{i} = (A+A')/2;
  toc
  [~,p]=chol(Matrices{i});
  if p~=0, error('not pos def. cond too high?'); end
  Bs{i} = randn(n,1);
end
%%
save matrices_n10000