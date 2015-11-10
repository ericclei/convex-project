function Atilde = LinearShrinkage(A, max_eig, min_eig, threshold, alpha)
    gamma = (1-alpha) * (max_eig - min_eig * threshold) / (alpha * (threshold - 1));
    n = size(A,1);
    Atilde = (1-alpha) * A + spdiags(ones(n,1)*alpha*gamma,0,n,n);
end

