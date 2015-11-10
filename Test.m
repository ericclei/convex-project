load('Matrices.mat')
threshold = .1/eps;
end_result = cell();
for i = 1:len(Matrices)
   A = Matrices{i};
   % Run SOR
   
   % Find obj value - x is result of linear system
   end_result{i,1} = x' * A * x - b' * x;
   
   % Find eigenvalues
   max_eig = eigs(A, 1);
   min_eig = eigs(A, 1, 'SM');
   shrinkage_result = zeros(10,1);
   for j = 1:10
       shA = LinearShrinkage(A, max_eig, min_eig, threshold, 0.01 * j);
       % Run SOR
       
       % Find obj value
       shrinkage_result(j) = x' * A * x - b' * x;
   end
   end_result{i,2} = shrinkage_result;
end
save('end_result.mat','end_result')