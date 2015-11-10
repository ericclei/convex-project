function Atilde = LinearShrinkage(A, max_eig, min_eig, threshold, alpha)
    gamma = (1-alpha) * (max_eig - min_eig * threshold) / (alpha * (threshold - 1));
    Atilde = (1-alpha) * A + alpha * gamma * eye(size(A));
end

