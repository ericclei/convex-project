clear
close all
clc
%%
load('matrices.mat')
threshold = .1/eps;
Er = cell(m,1);
EndResult = cell(m,1);
T = 10;
n_shrinkage = 50;
lo = .99;
hi = 1.01;
for k=1:m
  k
  A = Matrices{k};
  PerturbedMatrices = cell(T+1,1);
  end_result = cell(T+1,1+n_shrinkage);
  %%
  warning off
  for i = 1:T+1
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
     x = A\b;%x = successive_over_relaxation(A,b,w,1);
     toc
     % Find obj value - x is result of linear system
     end_result{i,1} = x;

     % Find eigenvalues
     max_eig = eigs(A, 1);
     min_eig = eigs(A, 1, 'SM');
     for j = 1:n_shrinkage
         shA = LinearShrinkage(A, max_eig, min_eig, threshold, 0.01 * j);
         % Run SOR
         x = shA\b;%x = successive_over_relaxation(shA,b,w,1);
         % Find obj value
         end_result{i,j+1}=x;
     end
  end
  warning on
  EndResult{k} = end_result;
  %%
  x_star = end_result{1,1};
  for i=1:T+1
    for j=1:1+n_shrinkage
      er(i,j) = norm(x_star-end_result{i,j});
    end
  end
  Er{k} = er;
end
%%
no_shrink_er = nan(m,1);
shrink_er = nan(m,n_shrinkage);
for k=1:m
  true_er = Er{k}(1,:);
  pert_er = mean(Er{k}(2:end,:));
  change = pert_er - true_er;
  no_shrink_er(k) = change(1);
  shrink_er(k,:) = change(2:end);
end
%%
figure,semilogy(.01:.01:.01*n_shrinkage,shrink_er')
xlabel('\alpha'),ylabel('L_2 error')
hold on
plot(.25, no_shrink_er,'x')
hold off
%%
save('midway_result.mat','Er')