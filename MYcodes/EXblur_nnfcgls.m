% EXblur_cgls_hybrid Example script, speckle deblurring problem

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

rng(0);  % Make sure this test is repeatable.

% Define the test problem.
NoiseLevel = 0.01;
options.trueImage = 'satellite';
[A, b, x, ProbInfo] = PRblurrotation(options);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);

% Run NN-FCGLS, use the true image to compute error norms, and find iteration
% where error is minimum (i.e., investigate semi-convergence).
options = IRset('x_true', x);
tic;
[X_nnfcgls, IterInfo_nnfcgls] = IRnnfcgls(A, bn, options);
time_nnfcgls = toc;

% Now use CGLS with the discrepancy principle as a stopping criterion.
% Use a large safety factor eta to simulate a situation where the noise
% level is quite uncertain.
options = IRset(options, 'RegParam', 'discrep', 'NoStop', 'off', 'NoiseLevel', NoiseLevel);
tic;
[X_nnfcgls_dp, IterInfo_nnfcgls_dp] = IRnnfcgls(A, bn, options);
time_nnfcgls_dp = toc;


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
    subplot(3,3,2), semilogy(1:200, IterInfo_nnfcgls.Enrm, 'b-', 'LineWidth', 1.5)
    hold on
    hl = legend('nnfcgls','location','northeast');
    set(hl,'interpreter','latex','fontsize',12)
    semilogy(IterInfo_nnfcgls.BestReg.It, IterInfo_nnfcgls.BestReg.Enrm, 'ro', 'LineWidth', 1.5, 'MarkerSize', 6)
    semilogy(IterInfo_nnfcgls_dp.its, IterInfo_nnfcgls_dp.Enrm(end), 'ms', 'LineWidth', 1.5, 'MarkerSize', 6)
    title('Error history','interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,5), imagesc(reshape(IterInfo_nnfcgls.BestReg.X, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(['Best NNFCGLS sol., $k$ = ' num2str(IterInfo_nnfcgls.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
    subplot(3,3,3), imagesc(reshape(X_nnfcgls_dp, ProbInfo.xSize)), axis off, axis image, colormap(gray)
    title(['DP NNFCGLS sol., $k$ = ',num2str(IterInfo_nnfcgls_dp.StopReg.It)],...
    'interpreter','latex','fontsize',18)
    set(gca,'fontsize',10)
    %
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
    semilogy(1:200, IterInfo_nnfcgls.Enrm, 'b-', 'LineWidth', LW)
    hold on
    semilogy(IterInfo_nnfcgls.BestReg.It, IterInfo_nnfcgls.BestReg.Enrm, 'ro', 'LineWidth', LW, 'MarkerSize', MS)
    semilogy(IterInfo_nnfcgls_dp.its, IterInfo_nnfcgls_dp.Enrm(end), 'ms', 'LineWidth', LW, 'MarkerSize', MS)
    hl = legend('{\tt IRnnfcgls} errors','optimal {\tt IRnnfcgls} stopping iteration', ...
        '{\tt IRnnfcgls} DP stopping iteration');
    set(hl,'interpreter','latex','fontsize',18)
    %
    figure(4), clf
    PRshowx(IterInfo_nnfcgls.BestReg.X, ProbInfo)
    title(['Best NNFCGLS sol., $k$ = ' num2str(IterInfo_nnfcgls.BestReg.It)],...
    'interpreter','latex','fontsize',18)
    %
    figure(5), clf
    PRshowx(X_nnfcgls_dp, ProbInfo)
    title(['DP NNFCGLS sol., $k$ = ',num2str(IterInfo_nnfcgls_dp.StopReg.It)],...
    'interpreter','latex','fontsize',18)
end


