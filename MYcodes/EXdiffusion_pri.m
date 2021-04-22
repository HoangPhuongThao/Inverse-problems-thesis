%Using PRI and Restricted PRI on 2D image deblurring problem

clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 1.5;  % Plot line width
MS = 10; % Size of markers on plots

rng(0);  % Make sure this test is repeatable.

% Define the test problem.
NoiseLevel = 0.01;
[A, b, x, ProbInfo] = PRdiffusion;
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
alpha = 1.01;
beta = 1e-4;
l = 15;
maxit = 100;
method = "rrgmres";
tol = 1e-2;

[X, rho_R, eta_X, outit, innit, time] = AFHpri(A,bn,alpha,NoiseLevel,l, maxit, method);
[X_tol, rho_R_tol, eta_X_tol, outit_tol, innit_tol, time_tol] = AFHpri(A,bn,alpha,NoiseLevel,l, maxit, method);

rel_err = zeros(1, outit);
nx = norm(x);
for i=1:outit
    rel_err(i) = norm(x-X(:,i))/nx;
end

rel_err_tol = zeros(1, outit_tol);
for i=1:outit_tol
    rel_err_tol(i) = norm(x-X_tol(:,i))/nx;
end

minRelErrInd = find(rel_err == min(rel_err));
minRelErrIndTol = find(rel_err_tol == min(rel_err_tol));

%if discrepancy principle fails to stop iterations soon enough, try finding
%the index of the first nondecreasing value in rho_R, rho_r2
tau = 1e-3;
for stopInd=1:length(rho_R)-1
    if (abs(rho_R(stopInd)- rho_R(stopInd+1))/rho_R(stopInd)) < tau
        break;
    end
    if stopInd==maxit-1
        stopInd=maxit;
    end
end
for stopIndTol=1:length(rho_R_tol)-1
    if (abs(rho_R_tol(stopIndTol)- rho_R_tol(stopIndTol+1))/rho_R_tol(stopIndTol)) < tau
        break;
    end
    if stopIndTol==maxit-1
        stopIndTol=maxit;
    end
end

figure(1), clf
PRshowx(x,ProbInfo)
title('True solution','interpreter','latex','fontsize',18)
set(gca,'fontsize',24)

figure(2), clf
PRshowb(bn,ProbInfo)
title('Noisy data','interpreter','latex','fontsize',18)
set(gca,'fontsize',24)

figure(3), clf
axes('FontSize', 24), hold on
semilogy(1:outit, rel_err, 'b-', 'LineWidth', LW)
hold on
semilogy(1:outit_tol, rel_err_tol, 'm-', 'LineWidth', LW)
semilogy(minRelErrInd, rel_err(minRelErrInd), 'ro', 'LineWidth', LW, 'MarkerSize', MS)
semilogy(minRelErrIndTol, rel_err_tol(minRelErrIndTol), 'ro', 'LineWidth', LW, 'MarkerSize', MS)
semilogy(stopInd, rel_err(stopInd), 'go', 'LineWidth', LW, 'MarkerSize', MS)
semilogy(stopIndTol, rel_err_tol(stopIndTol), 'gd', 'LineWidth', LW, 'MarkerSize', MS)
legend('PRI', 'PRI tol', 'minRelErr', 'minRelErrTol', 'stopped', 'stoppedTol')
xlabel('k iterations')
ylabel('||x_exact - x_k|| / ||x_exact||')
title('PRI Relative errors')

figure(4), clf
PRshowx(X(:,minRelErrInd), ProbInfo)
title(sprintf('Best Projected %s solution, k = %d iterations', method, minRelErrInd))

figure(5), clf
PRshowx(X_tol(:,minRelErrIndTol), ProbInfo)
title(sprintf('Best PRI tol %s solution, k = %d iterations', method, minRelErrIndTol))

figure(6), clf
PRshowx(X(:,stopInd), ProbInfo)
title(sprintf('Stopped PRI solution, k = %d iterations', stopInd))

figure(7), clf
PRshowx(X_tol(:,stopIndTol), ProbInfo)
title(sprintf('Stopped PRI tol solution, k = %d iterations', stopIndTol))