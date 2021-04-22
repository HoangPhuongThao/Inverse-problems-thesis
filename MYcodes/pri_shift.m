function [X, rho_R, eta_X, outit, innit, time] = pri_shift(A,b,alpha,delta,l, maxit, method)
%X - holds all approximate solutions generated each iteration
%delta - estimate of the norm of the noise
%alpha>=1 specified constant/safety factor
%l - maximum number of consecutive inner iterations
%maxit - maximum number of outer iterations (?)
%We'll start initiating from 1 instead of 0 as it was in the article
[m,n] = size(A);
k = 1;
X(:,k) = zeros(m,1);
R(:,k) = b;
dp = alpha*delta;
innit = 0;

tic;
%check the stopping criteion: 
while (norm(R(:,k)) > dp) && (k ~= maxit) 
    %reziduova rovnice - apply RRGMRES, CGLS or LSQR
    if method == "lsqr"
        [W,rho,eta] = lsqr_b(A,R(:,k),l);
    elseif method == "cgls"
        [W,rho,eta] = cgls(A,R(:,k),l);
    elseif method == "gmres"
        [W,rho,eta] = gmres(A,R(:,k),l);
    else
        [W,rho,eta] = rrgmres(A,R(:,k),l);
    end
    
    %The solution norm and residual norm are returned in eta and rho, resp.
    %discrepancy principle: ||r_k|| <= alpha*delta
    %terminate when the discr.pr. is satisfied or when nb it = l
    ix = find(rho <= dp, 1, 'first');
    if isempty(ix) == 1
        ix = l;
    end
    innit = innit + ix;
    
    % update the iterate
    x_new = X(:,k) + W(:,ix);
    % shift vector x_new by the min_val to the right --> x_new becomes
    % nonnegative
    min_val = min(x_new);
    if min_val < 0
        x_new = x_new - min_val*ones(size(x_new));
    end
    X(:,k+1) = x_new;
    k = k+1;
    R(:,k) = b- A*X(:,k);
end
time = toc;

outit = k;
rho_R = vecnorm(R(:, 1:outit));
eta_X = vecnorm(X(:, 1:outit));

end