function [ alpha, beta ] = Lanczos(A, w)
    n = length(A);
    v = zeros(n, 1);
    beta = zeros(n - 1, 1);
    alpha = zeros(n, 1);
    k = 0;
    while k == 0 || abs(beta(k)) > 1e-5 && k < n
        if k ~= 0
           temp = w;
           w = v / beta(k);
           v = - beta(k) * temp;
        end
        v = v + A * w;
        k = k + 1;
        alpha(k) = dot(w, v);
        v = v - alpha(k) * w;
        beta(k) = norm(v);
    end
    beta = beta(1:k - 1);
    alpha = alpha(1:k);
end

