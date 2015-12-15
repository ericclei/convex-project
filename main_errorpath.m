clear
close all
clc
%%
load('matrices_n10000.mat')
threshold = .01/eps;
Er = cell(m,1);
EndResult = cell(m,1);
T = 10;
n_shrinkage = 9;
n_methods = n_shrinkage+2;
ErrPath = cell(m,T+1,n_methods);
pert=.00001;
lo = 1-pert;
hi = 1+pert;
w = 1;
%%
for k=1:m
  k
  A = Matrices{k};
  x_star = A\b;
  original_A = A;
  PerturbedMatrices = cell(T,1);
  end_result = cell(T+1,n_methods);
  %%
  warning off
  er = nan(T+1,n_methods);
  for i = 1:T+1
    i
    if i>1
      A=original_A;
      P=(rand(nnz(A),1))*(hi-lo)+lo;
      nz = A~=0;
      NewA = A(nz).*P;
      A(nz) = NewA;
      PerturbedMatrices{i-1} = A;
    end
    % Run SOR
    tic
    [x,n_iter,ErrPath{k,i,1}] = symmetric_successive_over_relaxation(A,b,w,1,x_star);
    toc
    n_iter
    end_result{i,1} = x;
    er(i,1) = norm(x-x_star);
    % Run preconditioned SOR
    tic
    invP = spdiags(1./diag(A),0,n,n);
    [x,n_iter,ErrPath{k,i,2}] = symmetric_successive_over_relaxation(invP*A,invP*b,w,1,x_star);
    toc
    n_iter
    end_result{i,2} = x;
    er(i,2) = norm(x-x_star);
%%
    % Find eigenvalues
    max_eig = eigs(A, 1);
    min_eig = eigs(A, 1, 'SM');
    for j = 1:n_shrinkage
       j
       A_ls = LinearShrinkage(A, max_eig, min_eig, threshold, 0.1 * j);
       % Run SOR
       tic
       [x,n_iter,ErrPath{k,i,j+n_methods-n_shrinkage}] = symmetric_successive_over_relaxation(A_ls,b,w,1,x_star);
       toc
       n_iter
       % Find obj value
       end_result{i,j+n_methods-n_shrinkage}=x;
       er(i,j+n_methods-n_shrinkage) = norm(x-x_star);
    end
  end
  warning on
  EndResult{k} = end_result;
  Er{k} = er;
end
%%
no_shrink_er = nan(m,1);
precond_er = nan(m,1);
shrink_er = nan(m,n_shrinkage);
for k=1:m
  true_er = Er{k}(1,:);
  pert_er = mean(Er{k}(2:end,:),1);
  change = pert_er - true_er;
  no_shrink_er(k) = change(1);
  precond_er(k) = change(2);
  shrink_er(k,:) = change(1+n_methods-n_shrinkage:end);
end
%%
figure,semilogy(.01:.01:.01*n_shrinkage,shrink_er')
xlabel('\alpha'),ylabel('L_2 error')
hold on
plot(.25, no_shrink_er,'x')
hold off
%%
clear Matrices
save('midway_errorpath.mat')
% Er: mx1 cell
% Er{k} is a (T+1) x (n_shrinkage+2) matrix of error values (norm(Ax-b))
% - first row is original A, each next row is a different perturbation of A
% - first col is direct solver, second is preconditioned, rest are shrinkage with different parameters
% EndResult is the same as Er, except it contains cells of the solutions x instead of error
% ErrPath is the same, except it contains the error path of the solutions over iterations