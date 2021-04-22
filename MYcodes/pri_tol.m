function [X, rho_R, eta_X, outit, innit, time] = pri_tol(A,b,alpha,delta,l, maxit, method, tol, tau)
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
    for i=1:l-1
        if (abs(rho(i) - rho(i+1))/rho(i) < tau)
            break;
        end
    end
    innit = innit + min(ix,i);
    
    % update the iterate
    x_new = X(:,k) + W(:,ix);
    % project element to zero if it's smaller than -tol
    for i=1:m
        if x_new(i) < -tol
            x_new(i) = 0;
        end
    end
    X(:,k+1) = x_new;
    R(:,k+1) = b- A*X(:,k+1);
    if (abs(norm(R(:,k)) - norm(R(:,k+1)))/norm(R(:,k)) < tau)
        break;
    end
    k = k+1;
end
time = toc;

outit = k;
rho_R = vecnorm(R(:, 1:outit));
eta_X = vecnorm(X(:, 1:outit));

end