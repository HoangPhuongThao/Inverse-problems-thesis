%new version of PRI function for A function handle

function [X, rho_R, eta_X, outit, innit, time] = AFHpritol(A,b,alpha,delta,l, maxit, method, tol)
%X - holds all approximate solutions generated each iteration
%delta - estimate of the norm of the noise
%alpha>=1 specified constant/safety factor
%l - maximum number of consecutive inner iterations
%maxit - maximum number of outer iterations (?)
%We'll start initiating from 1 instead of 0 as it was in the article
m = length(b);
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
        [W,info] = IRhybrid_lsqr(A,R(:,k),1:l);
    elseif method == "cgls"
        [W,info] = IRcgls(A,R(:,k),1:l);
    else
        [W,info] = IRrrgmres(A,R(:,k),1:l);
    end
    
    %The solution norm and residual norm are returned in eta and rho, resp.
    %discrepancy principle: ||r_k|| <= alpha*delta
    %terminate when the discr.pr. is satisfied or when nb it = l
    ix = find(info.Rnrm <= dp, 1, 'first');
    if isempty(ix) == 1
        ix = l;
    end
    innit = innit + ix;
    
    % update the iterate
    x_new = X(:,k) + W(:,ix);
    % project element to zero if it's smaller than -tol
    for i=1:m
        if x_new(i) < -tol
            x_new(i) = 0;
        end
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