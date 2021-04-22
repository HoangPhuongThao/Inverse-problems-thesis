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
options.trueImage = 'satellite';
[A, b, x, ProbInfo] = PRblurrotation(options);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
alpha = 2;
beta = 1e-4;
l = 10;
maxit = 30;
method = "cgls";
tol = 0.001;
tau = 1e-3;

[X, rho_R, eta_X, outit, innit, time] = projected_rest_it(A,bn,alpha,NoiseLevel,l, maxit, method, tau);
% [X, rho_R, eta_X, outit, innit, time] = pri_shift(A,bn,alpha,NoiseLevel,l, maxit, method);
[X_tol, rho_R_tol, eta_X_tol, outit_tol, innit_tol, time_tol] = pri_tol(A,bn,alpha,NoiseLevel,l, maxit, method, tol, tau);
[X2, rho_R2, eta_X2, outit2, innit2, time2, m] = restrictedPRI(A,bn,alpha,NoiseLevel,l, maxit,beta, method, tau);

rel_err = zeros(1, outit);
nx = norm(x);
for i=1:outit
    rel_err(i) = norm(x-X(:,i))/nx;
end

rel_err_tol = zeros(1, outit_tol);
for i=1:outit_tol
    rel_err_tol(i) = norm(x-X_tol(:,i))/nx;
end

rel_err2 = zeros(1, outit2);
for i=1:outit2
    rel_err2(i) = norm(x-X2(:,i))/nx;
end

minRelErrInd = find(rel_err == min(rel_err));
minRelErrIndTol = find(rel_err_tol == min(rel_err_tol));
minRelErrInd2 = find(rel_err2 == min(rel_err2));

rho_R = rho_R/norm(bn);
rho_R2 = rho_R2/norm(bn);
rho_R_tol = rho_R_tol/norm(bn);

%if discrepancy principle fails to stop iterations soon enough, try finding
%the index of the first nondecreasing value in rho_R, rho_r2
% tau = 1e-3;
% for stopInd=1:length(rho_R)-1
%     if (abs(rho_R(stopInd)- rho_R(stopInd+1))/rho_R(stopInd)) < tau
%         break;
%     end
%     if stopInd==maxit-1
%         stopInd=maxit;
%     end
% end
% for stopInd2=1:length(rho_R2)-1
%     if (abs(rho_R2(stopInd2)- rho_R2(stopInd2+1))/rho_R2(stopInd2)) < tau
%         break;
%     end
%     if stopInd2==maxit-1
%         stopInd2=maxit;
%     end
% end
% for stopIndTol=1:length(rho_R_tol)-1
%     if (abs(rho_R_tol(stopIndTol)- rho_R_tol(stopIndTol+1))/rho_R_tol(stopIndTol)) < tau
%         break;
%     end
%     if stopIndTol==maxit-1
%         stopIndTol=maxit;
%     end
% end

if strcmp(dispres, 'subplots')
    subplot(2,3,1), imagesc(reshape(x, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title('True solution','interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(2,3,2), imagesc(reshape(bn, ProbInfo.bSize)), axis off, axis image, colormap(gray)
    title('Noisy data','interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(2,3,3)
    imagesc(reshape(X(:,minRelErrInd), ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(sprintf('Best PRI solution, k = %d iterations', minRelErrInd))
    set(gca,'fontsize',10)
    %
    subplot(2,3,4)
    imagesc(reshape(X(:,stopInd), ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(sprintf('Stopped PRI solution, k = %d iterations', stopInd))
    set(gca,'fontsize',10)
    %
    subplot(2,3,5)
    imagesc(reshape(X2(:,minRelErrInd2), ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(sprintf('Best RPRI solution, k = %d iterations', minRelErrInd2))
    set(gca,'fontsize',10)
    %
    subplot(2,3,6)
    imagesc(reshape(X2(:,stopInd2), ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(sprintf('Stopped RPRI solution, k = %d iterations', stopInd2))
    set(gca,'fontsize',10)
    %
    figure
    semilogy(1:outit, rel_err, 'b-', 'LineWidth', LW)
    hold on
    semilogy(1:outit2, rel_err2, 'r-', 'LineWidth', LW)
    semilogy(minRelErrInd, rel_err(minRelErrInd), 'bo', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(minRelErrInd2, rel_err2(minRelErrInd2), 'rx', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(stopInd, rel_err(stopInd), 'go', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(stopInd2, rel_err2(stopInd2), 'gx', 'LineWidth', LW, 'MarkerSize', MS)
    legend('PRI','restricted PRI', 'minRelErr', 'minRelErr2', 'stopped', 'stopped2')
    xlabel('k iterations')
    ylabel('||x_exact - x_k|| / ||x_exact||')
    title('Relative error')
    set(gca,'fontsize',10)
elseif strcmp(dispres, 'manyplots')
    figure(1), clf
    PRshowx(x, ProbInfo)
    set(gca,'fontsize',24)
    title('True solution','interpreter','latex','fontsize',18)

    figure(2), clf
    PRshowb(b, ProbInfo)
    set(gca,'fontsize',24)
    title('Noisy data','interpreter','latex','fontsize',18)

    figure(3), clf
    axes('FontSize', 20), hold on
    semilogy(1:outit, rel_err, 'b-', 'LineWidth', LW)
    hold on
    semilogy(1:outit2, rel_err2, 'r-', 'LineWidth', LW)
    semilogy(1:outit_tol, rel_err_tol, 'm-', 'LineWidth', LW)
    semilogy(minRelErrInd, rel_err(minRelErrInd), 'bo', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(minRelErrInd2, rel_err2(minRelErrInd2), 'rx', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(minRelErrIndTol, rel_err_tol(minRelErrIndTol), 'md', 'LineWidth', LW, 'MarkerSize', MS)
%     semilogy(stopInd, rel_err(stopInd), 'go', 'Line   Width', LW, 'MarkerSize', MS)
%     semilogy(stopInd2, rel_err2(stopInd2), 'gx', 'LineWidth', LW, 'MarkerSize', MS)
%     semilogy(stopIndTol, rel_err_tol(stopIndTol), 'gd', 'LineWidth', LW, 'MarkerSize', MS)
%     legend('PRI','restricted PRI', 'PRI tol', 'minRelErr', 'minRelErr2', 'minRelErrTol', 'stopped', 'stopped2', 'stoppedTol')
    legend('PRI','RSPRI', 'MPRI', 'minRelErrPRI', 'minRelErrRSPRI', 'minRelErrMPRI')
    xlabel('k')
    ylabel('||x^{exact} - x_k|| / ||x^{exact}||')
    title('Relative error')

    figure(4), clf
    PRshowx(X(:,minRelErrInd), ProbInfo)
    title(sprintf('PRI %s solution, k = %d iterations', method, minRelErrInd))

    figure(5), clf
    PRshowx(X2(:,minRelErrInd2), ProbInfo)
    title(sprintf('RSPRI %s solution, k = %d iterations', method, minRelErrInd2))
    
    figure(6), clf
    PRshowx(X_tol(:,minRelErrIndTol), ProbInfo)
    title(sprintf('MPRI %s solution, k = %d iterations', method, minRelErrIndTol))
    
%     figure(7), clf
%     PRshowx(X(:,stopInd), ProbInfo)
%     title(sprintf('Stopped PRI solution, k = %d iterations', stopInd))
%     
%     figure(8), clf
%     PRshowx(X2(:,stopInd2), ProbInfo)
%     title(sprintf('Stopped restricted PRI solution, k = %d iterations', stopInd2))
%     
%     figure(9), clf
%     PRshowx(X_tol(:,stopIndTol), ProbInfo)
%     title(sprintf('Stopped PRI tol solution, k = %d iterations', stopIndTol))
end
