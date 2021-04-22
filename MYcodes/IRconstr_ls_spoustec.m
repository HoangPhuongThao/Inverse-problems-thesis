clear
clc

LW = 2;  % Plot line width.
MS = 10; % Size of markers on plots.

% 1D problem - shaw
% [A,b,x] = shaw(100);
% e = randn(size(b)); e = e/norm(e); e = 1e-4*norm(b)*e;
% bn = b + e;

% % Define test problem - diffusion
% n = 64;                                   % Problem size.
% NoiseLevel = 0.005;                       % Relative noise level in data.
% [A,b,x,ProbInfo] = PRdiffusion(n);        % Get the test problem.
% [bn,NoiseInfo] = PRnoise(b, NoiseLevel); 

% EXblur problem
NoiseLevel = 0.01;
options.trueImage = 'dotk';
[A, b, x, ProbInfo] = PRblurspeckle(options);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);

% Run CGLS, use the true image to compute error norms, and find iteration
% where error is minimum (i.e., investigate semi-convergence).
options = IRset('x_true', x);

% Call the Projected restarted iteration
options.MaxIterIn = 15;
options.MaxIterOut = 100;
options.RegParam = 'discrep';
options.NoiseLevel = NoiseLevel;
tic;
[X, Info] = IRconstr_ls(A, b, options);
time = toc;

% Approximate solution
figure(1)
PRshowx(Info.StopReg.X, ProbInfo);
title(['IRconstr_ls solution for k = ', num2str(Info.its)],'interpreter','latex','fontsize',20)
set(gca,'fontsize',13)

%Approximate solution for 1D problem
% figure(1)
% hold on
% plot(x, 'b'); plot(X, 'g');
% legend('exact solution', 'pri')

% Relative residual history
figure(2)
res = zeros(1,size(Info.itsInOut,1));
for i = 1:size(Info.itsInOut,1)
    tot_it = Info.itsInOut(i,3);
    res(i) = Info.Rnrm(tot_it);
end
semilogy(res,'-','linewidth',LW)
title('Residual history','interpreter','latex','fontsize',22);
set(gca,'fontsize',13)
xlabel('k iterations')
ylabel('||b-Ax_k||/||b||')

% Relative error history
figure(3)
normx = norm(x);
Enrm = zeros(1,size(Info.Xout,2));
for i = 1:size(Info.Xout,2)
    Enrm(i) = norm(x - Info.Xout(:,i))/ normx;
end
semilogy(Enrm,'-','linewidth',LW)
title('Error history','interpreter','latex','fontsize',22);
set(gca,'fontsize',13)
xlabel('k iterations')
ylabel('||x_{exact} - x_k|| / ||x_{exact}||')
