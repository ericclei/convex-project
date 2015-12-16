function [x,err_path] = conjugate_gradient_shrinkage(A,b,eta,original_A)
    d = size(A,1);
    x = zeros(d,1);
    r=b-A*x;
    if nargout > 1
      err_path = nan(d+1,1);
      err_path(1) = norm(b-original_A*x);
    end
    tol = 1e-10;
    p=r;
    rsold=r'*r;
    if sqrt(rsold/d)<tol
      return
    end
    num_iter = 0;
    
    for i=1:length(b)
        num_iter = num_iter+1;
        Ap=A*p;
        old_denom = p'*Ap;
        denom = old_denom + eta*sum(p.^2);
        alpha=rsold/denom;
        x=x+alpha*p;
        if nargout > 1
          err_path(1+num_iter) = norm(b-original_A*x);
        end
        r=r-alpha*Ap;
        rsnew=r'*r;
        if mod(i,1000)==0
          disp(i)
          disp(sqrt(rsnew/d))
        end
        if sqrt(rsnew/d)<tol
              break;
        end
        p=r+(rsnew/rsold)*p;
        rsold=rsnew;
    end
end