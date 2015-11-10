function [ min_lambda, max_lambda ] = FindEigenValues( A )
    w = rand(length(A), 1);
    w = w / norm(w);
    [a, b] = Lanczos(A, w);
%     e = trideig(a, b)
    A_trid = diag(a, 0) + diag(b, -1) + diag(b, 1);
    e = eig(A_trid);
%     e = trideigs(a, b)
    min_lambda = min(e);
    max_lambda = max(e);
end

