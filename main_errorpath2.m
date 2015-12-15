clear
close all
clc
%%
load('matrices.mat')
threshold = .1/eps;
Er = cell(m,1);
EndResult = cell(m,1);
T = 1;
n_shrinkage = 1;
n_methods = n_shrinkage+2;
ErrPath = cell(m,T+1,n_methods);
pert = .0001;
lo = 1-pert;
hi = 1+pert;
w = .7;
z=2;
for k=z
  k
  A = Matrices{k};
  x_star = A\b;
  PerturbedMatrices = cell(T+1,1);
  end_result = cell(T+1,n_methods);
  %%
  warning off
  er = nan(T+1,n_methods);
  for i = 1:T+1
     i
     if i>1
       A=PerturbedMatrices{1};
       P=(rand(nnz(A),1)-lo)/(hi-lo);
       nz = A~=0;
       NewA = A(nz).*P;
       A(nz) = NewA;
     end
     PerturbedMatrices{i} = A;
     % Run SOR
     tic
     [x,~,ErrPath{k,i,1}] = symmetric_successive_over_relaxation(A,b,w,1,x_star);
%      x = jacobi_method(A,b);
     toc
     end_result{i,1} = x;
%      x_star = end_result{1,1};
     er(i,1) = norm(x-x_star);
     % Run preconditioned SOR
     tic
     invP = spdiags(1./diag(A),0,n,n);
     [x,~,ErrPath{k,i,2}] = symmetric_successive_over_relaxation(invP*A,invP*b,w,1,x_star);
%      x = jacobi_method(invP*A,invP*b);
     toc
     end_result{i,2} = x;
     er(i,2) = norm(x-x_star);
%%
     % Find eigenvalues
     max_eig = eigs(A, 1);
     min_eig = eigs(A, 1, 'SM');
     for j = 1:n_shrinkage
         j
         shA = LinearShrinkage(A, max_eig, min_eig, threshold, 0.1 * j);
         % Run SOR
         [x,~,ErrPath{k,i,j+n_methods-n_shrinkage}] = symmetric_successive_over_relaxation(shA,b,w,1,x_star);
%          x = jacobi_method(shA,b);
         % Find obj value
         end_result{i,j+n_methods-n_shrinkage}=x;
         er(i,j+n_methods-n_shrinkage) = norm(x-x_star);
     end
  end
  warning on
  EndResult{k} = end_result;
  Er{k} = er;
  %%
%   x_star = end_result{1,1};
%   for i=1:T+1
%     for j=1:1+n_shrinkage
%       er(i,j) = norm(x_star-end_result{i,j});
%     end
%   end
end
%%
no_shrink_er = nan(m,1);
shrink_er = nan(m,n_shrinkage);
for k=z
  true_er = Er{k}(1,:);
  pert_er = mean(Er{k}(2:end,:));
  change = pert_er - true_er;
  no_shrink_er(k) = change(1);
  shrink_er(k,:) = change(2:end);
end
%%
save(sprintf('midway_errorpath%d.mat',z),'Er','EndResult','ErrPath')