function [x] = conjugate_gradient(A,b,x0)
    d = size(A,1);
    if ~exist('x0','var')
      x0 = zeros(d,1);
    end
    x = x0;
    r=b-A*x;
    p=r;
    rsold=r'*r;
    if sqrt(rsold/d)<1e-10
      return
    end

    for i=1:length(b)
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        disp(sqrt(rsnew/d))
        if sqrt(rsnew/d)<1e-10
              break;
        end
        p=r+(rsnew/rsold)*p;
        rsold=rsnew;
    end
end