function D = threshold_diagonal(D,L)
% put L as min value of diagonal matrix D
D = diag(max(L,diag(D)));