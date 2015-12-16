clear
close all
clc
%%
load('matrices_n10000.mat')
Er = cell(m,1);
EndResult = cell(m,1);
T = 1;
pert=.00001;
lo = 1-pert;
hi = 1+pert;
etas = 10^-5;
n_shrinkage = length(etas);
n_methods = n_shrinkage+1;
PerturbedMatrices = cell(m,T);
%%
for k=1:m
  k
  A = Matrices{k};
  original_A = A;
  b = Bs{k};
  x_star = A\b;
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
       PerturbedMatrices{k,i-1} = A;
     end
     % Exact solution
     tic
     x = A\b;
     toc
     end_result{i,1} = x;
     x_star = end_result{1,1};
     er(i,1) = norm(original_A*x-b);%norm(x-x_star);
     % Conj grad 
%      tic
     [x,err_path_plain] = conjugate_gradient(A,b,original_A);
%      toc

     % Conj grad shrinkage
     for j = 1:n_shrinkage
         j
         tic
         [x,err_path] = conjugate_gradient_shrinkage(A,b,etas(j),original_A);
         toc
         end_result{i,j+1}=x;
         er(i,j+1) = norm(original_A*x-b);%norm(x-x_star);
     end
  end
  warning on
  EndResult{k} = end_result;
  Er{k} = er;
end
%%
figure,semilogy([err_path_plain err_path],'linewidth',3)
xlim([0 1e4])
xlabel('Iteration')
ylabel('Error')
legend({'conjugate gradient' 'spectral filtering conjugate gradient'})
set(gca,'FontSize',20)
%%
no_shrink_er = nan(m,1);
precond_er = nan(m,1);
shrink_er = nan(m,n_shrinkage);
OrigEr=nan(m,n_methods);
PertEr=nan(m,T,n_methods);
for k=1:m
  orig_er = Er{k}(1,:);
  pert_er = Er{k}(2:end,:);
  mean_pert_er = mean(pert_er);
  OrigEr(k,:) = orig_er;
  PertEr(k,:,:) = pert_er;
  change = mean_pert_er - orig_er;
  no_shrink_er(k) = change(1);
  precond_er(k) = change(2);
  shrink_er(k,:) = change(1+n_methods-n_shrinkage:end);
end
% %%
% BestEtaIdx=nan(m,T);
% for k=1:m
%   for i=1:T
%     [~,BestEtaIdx(k,i)] = min(PertEr(k,i,1+(n_methods-n_shrinkage):end));
%   end
% end
%%
figure
hold on
% for k=1:m
%   if k==3, continue; end
%   semilogx(mean(PertEr(k,:,1)),mean(PertEr(k,:,47)),'o') %eta=10^-5
% end
idx = [1:2 4:m];
scatter(squeeze(PertEr(idx,1,1)),squeeze(PertEr(idx,:,2)),'filled');
% plot([.75 1.25],[.75 1.25],'col','black')
hold off
xlabel('L_2 error plain')
ylabel('L_2 error with SF')
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
%%
% no_shrink_er = nan(m,1);
% shrink_er = nan(m,n_shrinkage);
% for k=1:m
%   true_er = Er{k}(1,:);
%   pert_er = mean(Er{k}(2:end,:));
%   change = pert_er - true_er;
%   no_shrink_er(k) = change(1);
%   shrink_er(k,:) = change(2:end);
% end
% %%
% figure,semilogy(.01:.01:.01*n_shrinkage,shrink_er')
% xlabel('\alpha'),ylabel('L_2 error')
% hold on
% plot(.25, no_shrink_er,'x')
% hold off
%%
save('conjgrad_result.mat');