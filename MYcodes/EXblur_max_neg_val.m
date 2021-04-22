% EXblur_cgls_hybrid Example script, speckle deblurring problem

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
dispres = 'subplots';
% dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

rng(0);  % Make sure this test is repeatable.

% Define the test problem.
NoiseLevel = 0.01;
options.trueImage = 'satellite';
[A, b, x, ProbInfo] = PRblurrotation(options);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);

% Run CGLS, use the true image to compute error norms, and find iteration
% where error is minimum (i.e., investigate semi-convergence).
options = IRset('x_true', x);
tic;
[X, X_all, IterInfo_cgls] = IRcgls2(A, bn, options);
time_cgls = toc;

% Now use CGLS with the discrepancy principle as a stopping criterion.
% Use a large safety factor eta to simulate a situation where the noise
% level is quite uncertain.
options = IRset(options, 'NoiseLevel', NoiseLevel);
tic;
[X_cgls_dp, IterInfo_cgls_dp] = IRcgls(A, bn, options);
time_cgls_dp = toc;

% Run the Hybrid LSQR method, with NoStop on.
options = IRset(options, 'RegParam', 'wgcv', 'NoStop', 'on');
% options = IRset(options, 'NoiseLevel', NoiseLevel);
tic;
[X_hybrid, X_all_lsqr, IterInfo_hybrid] = IRhybrid_lsqr2(A, bn, options);
time_hybrid = toc;

% Now run RRGMRES.
eta = 1.01;  % Safety factor.
K = 1:100;   % Iterations.
options = IRset('x_true', x, 'NoStop', 'on', 'NoiseLevel', NoiseLevel, 'eta', eta);
[X_rrgmres, IterInfo_rrgmres] = IRrrgmres(A,bn,K,options);


% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends
if strcmp(dispres, 'subplots')
    figure(1), clf
    subplot(3,3,1), imagesc(reshape(x, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title('True solution','interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,4), imagesc(reshape(bn, ProbInfo.bSize)), axis off, axis image, colormap(gray)
    title('Noisy data','interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,2), semilogy(1:100, IterInfo_cgls.Enrm, 'b-', 'LineWidth', 1.5)
    hold on
    semilogy(0:100, [norm(bn); IterInfo_hybrid.Enrm], 'k-.', 'LineWidth', 1.5)
    hl = legend('cgls','hybrid','location','northeast');
    set(hl,'interpreter','latex','fontsize',12)
    semilogy(IterInfo_cgls.BestReg.It, IterInfo_cgls.BestReg.Enrm, 'ro', 'LineWidth', 1.5, 'MarkerSize', 6)
    semilogy(IterInfo_cgls_dp.its, IterInfo_cgls_dp.Enrm(end), 'ms', 'LineWidth', 1.5, 'MarkerSize', 6)
    semilogy(IterInfo_hybrid.StopReg.It, IterInfo_hybrid.Enrm(IterInfo_hybrid.StopReg.It), 'mx', 'LineWidth', LW, 'MarkerSize', MS)
    title('Error history','interpreter','latex','fontsize',18)
    axis([0,100,1.5e-1,IterInfo_hybrid.Enrm(1)])
    set(gca,'fontsize',10)
    %
    subplot(3,3,5), imagesc(reshape(IterInfo_cgls.BestReg.X, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(['Best CGLS sol., $k$ = ' num2str(IterInfo_cgls.BestReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,3), imagesc(reshape(X_cgls_dp, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(['DP CGLS sol., $k$ = ',num2str(IterInfo_cgls_dp.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,6), imagesc(reshape(IterInfo_hybrid.StopReg.X, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(['DP IRhybrid\_lsqr sol., $k$ = ',num2str(IterInfo_hybrid.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
elseif strcmp(dispres, 'manyplots')
    figure(1), clf
    PRshowx(x, ProbInfo)
    set(gca,'fontsize',24)
    title('True solution','interpreter','latex','fontsize',18)
    %
    figure(2), clf
    PRshowb(b, ProbInfo)
    set(gca,'fontsize',24)
    title('Noisy data','interpreter','latex','fontsize',18)
    %
    figure(3), clf
    axes('FontSize', 24), hold on
    semilogy(1:100, IterInfo_cgls.Enrm, 'b-', 'LineWidth', LW)
    hold on
    semilogy(0:100, [norm(bn); IterInfo_hybrid.Enrm], 'k-.', 'LineWidth', LW)
    semilogy(IterInfo_cgls.BestReg.It, IterInfo_cgls.BestReg.Enrm, 'ro', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(IterInfo_cgls_dp.its, IterInfo_cgls_dp.Enrm(end), 'ms', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(IterInfo_hybrid.StopReg.It, IterInfo_hybrid.StopReg.Enrm, 'mx', 'LineWidth', LW, 'MarkerSize', MS)
    hl = legend('{\tt IRcgls} errors','{\tt IRhybrid\_lsqr} errors', ...
      'optimal {\tt IRcgls} stopping iteration','{\tt IRcgls} DP stopping iteration', ...
      '{\tt IRhybrid\_lsqr} DP stopping iteration');
    set(hl,'interpreter','latex','fontsize',18)
    title('Error history','interpreter','latex','fontsize',18)
    axis([0,100,0.15,IterInfo_hybrid.Enrm(1)])
    %
    figure(4), clf
    PRshowx(IterInfo_cgls.BestReg.X, ProbInfo)
    title(['Best CGLS sol., $k$ = ' num2str(IterInfo_cgls.BestReg.It)],...
    'interpreter','latex','fontsize',18)
    %
    figure(5), clf
    PRshowx(X_cgls_dp, ProbInfo)
    title(['DP CGLS sol., $k$ = ',num2str(IterInfo_cgls_dp.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    %
    figure(6)
    PRshowx(IterInfo_hybrid.StopReg.X, ProbInfo)
    title(['DP IRhybrid\_lsqr sol., $k$ = ',num2str(IterInfo_hybrid.StopReg.It)],...
    'interpreter','latex','fontsize',18)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display largest negative values at each iteration compared with relative
% errors (for CGLS, hybrid LSQR, RRGMRES)
figure, clf
axes('FontSize', 15), hold on
yyaxis right
ylabel('relative error')
semilogy(1:100, IterInfo_cgls.Enrm, 'r-', 'LineWidth', LW)
semilogy(IterInfo_cgls.BestReg.It, IterInfo_cgls.BestReg.Enrm, 'rx', 'LineWidth', LW, 'MarkerSize', MS)
hold on

cgls_largest_neg = zeros(size(X_all,2),1);
for k=1:size(X_all,2)
cgls_largest_neg(k) = min(X_all(:,k));
end
cgls_min = min(cgls_largest_neg);

yyaxis left
ylabel('largest negative value')
semilogy(1:100, cgls_largest_neg, 'b-', 'LineWidth', LW)
semilogy(find(cgls_largest_neg == cgls_min), cgls_min, 'bx', 'LineWidth', LW, 'MarkerSize', MS)
hl = legend('{\tt CGLS} largest neg. values', '{\tt CGLS} max neg. value',...
    '{\tt CGLS} relative errors','{\tt CGLS} min relative error');
set(hl,'interpreter','latex','fontsize',18)
xlabel('iterations')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, clf
axes('FontSize', 15), hold on
yyaxis right
semilogy(1:100, IterInfo_hybrid.Enrm, 'b-', 'LineWidth', LW)
semilogy(IterInfo_hybrid.BestReg.It, IterInfo_hybrid.BestReg.Enrm, 'bx', 'LineWidth', LW, 'MarkerSize', MS)
hold on

lsqr_largest_neg = zeros(size(X_all_lsqr,2),1);
for k=1:size(X_all_lsqr,2)
lsqr_largest_neg(k) = min(X_all_lsqr(:,k));
end
lsqr_min = min(lsqr_largest_neg);

yyaxis left
semilogy(1:100, lsqr_largest_neg, 'r-', 'LineWidth', LW)
semilogy(find(lsqr_largest_neg == lsqr_min), lsqr_min, 'rx', 'LineWidth', LW, 'MarkerSize', MS)
hl = legend('{\tt LSQR} relative errors','{\tt LSQR} min relative error',...
    '{\tt LSQR} largest neg. values', '{\tt LSQR} max neg. value');
set(hl,'interpreter','latex','fontsize',18)
xlabel('iterations')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, clf
axes('FontSize', 15), hold on
yyaxis right
semilogy(1:100, IterInfo_rrgmres.Enrm, 'b-', 'LineWidth', LW)
semilogy(IterInfo_rrgmres.BestReg.It, IterInfo_rrgmres.BestReg.Enrm, 'bx', 'LineWidth', LW, 'MarkerSize', MS)
hold on

rrgmres_largest_neg = zeros(size(X_rrgmres,2),1);
for k=1:size(X_rrgmres,2)
rrgmres_largest_neg(k) = min(X_rrgmres(:,k));
end
rrgmres_min = min(rrgmres_largest_neg);

yyaxis left
semilogy(1:100, rrgmres_largest_neg, 'r-', 'LineWidth', LW)
semilogy(find(rrgmres_largest_neg == rrgmres_min), rrgmres_min, 'rx', 'LineWidth', LW, 'MarkerSize', MS)
hl = legend('{\tt RRGMRES} relative errors','{\tt RRGMRES} min relative error',...
    '{\tt RRGMRES} largest neg. values', '{\tt RRGMRES} max neg. value');
set(hl,'interpreter','latex','fontsize',18)
xlabel('iterations')

return

% A number of instructions useful to save the displayed figures follow;
% the default is not to execute them. If you wish to save the displayed
% figures in the dedicated 'Results' folder, please comment the above
% return statement
oldcd = cd;
if strcmp(dispres, 'subplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -dpng -r300 EXblur_cgls_hybrid
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), print -depsc -r300 EXblur_cgls_hybrid_a
    figure(2), print -depsc -r300 EXblur_cgls_hybrid_b
    figure(3), print -depsc -r300 EXblur_cgls_hybrid_c
    figure(4), print -depsc -r300 EXblur_cgls_hybrid_d
    figure(5), print -depsc -r300 EXblur_cgls_hybrid_e
    figure(6), print -depsc -r300 EXblur_cgls_hybrid_f
end
cd(oldcd)

% Uncomment the following return statement if you wish to save the
% displayed figures as MATLAB figures

% return

oldcd = cd;
if strcmp(dispres, 'subplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    figure(1), saveas('EXblur_cgls_hybrid.fig')
elseif strcmp(dispres, 'manyplots')
    try
        cd('Results')
    catch
        mkdir('Results')
        cd('Results')
    end
    saveas(figure(1), 'EXblur_cgls_hybrid_a.fig')
    saveas(figure(2), 'EXblur_cgls_hybrid_b.fig')
    saveas(figure(3), 'EXblur_cgls_hybrid_c.fig')
    saveas(figure(4), 'EXblur_cgls_hybrid_d.fig')
    saveas(figure(5), 'EXblur_cgls_hybrid_e.fig')
    saveas(figure(6), 'EXblur_cgls_hybrid_f.fig')
end
cd(oldcd)