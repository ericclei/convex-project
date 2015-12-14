clear
close all
clc
%%
load('matrices.mat')
Er = cell(m,1);
EndResult = cell(m,1);
T = 1;
lo = .99;
hi = 1.01;
etas = 10.^(-50:0);
n_shrinkage = length(etas);
%%
for k=1:m
  k
  A = Matrices{k};
  PerturbedMatrices = cell(T+1,1);
  end_result = cell(T+1,1+n_shrinkage);
  %%
  warning off
  er = nan(T+1,1+n_shrinkage);
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
     % Exact solution
     tic
     x = A\b;%x = successive_over_relaxation(A,b,w,1);
     toc
     end_result{i,1} = x;
     x_star = end_result{1,1};
     er(i,1) = norm(x-x_star);
     % Conj grad 
%      tic
%      x = conjugate_gradient(A,b);
%      toc
%      end_result{i,2} = x;

     % Conj grad shrinkage
     for j = 1:n_shrinkage
         j
         tic
         x = conjugate_gradient_shrinkage(A,b,etas(j));
         toc
         end_result{i,j+1}=x;
         er(i,j+1) = norm(x-x_star)
     end
  end
  warning on
  EndResult{k} = end_result;
  Er{k} = er;
  %%
%   x_star = end_result{1,1};
%   er = nan(T+1, 1+n_shrinkage);
%   for i=1:T+1
%     for j=1:1+n_shrinkage
%       er(i,j) = norm(x_star-end_result{i,j});
%     end
%   end
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
save('conjgrad_result.mat','Er')