clear
clc

%This script runs tests on image deblurring problems.
%We compare the results of these 6 algorithms: CGLS, Projected CGLS,
%NNFCGLS, PRI, PRI_tol, RSPRI

true_img = {'dotk', 'satellite', 'hst'};
NoiseLevels = [1e-2; 1e-6];
alpha = [1.01; 1.1]; %safety factor
beta = 1e-4; %for Armijo's rule
l = 10; %max inner iter.
maxit = 30; %max outer iter.
method = "cgls";
tol = [1e-2; 1e-3; 1e-6]; %for PRI_tol
tau = [1e-2; 1e-3; 1e-6]; %for the stopping criteria

combinations = length(true_img)*length(NoiseLevels)*length(alpha)*length(tol)*length(tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_pri = zeros(combinations,1);
Res_pri = zeros(combinations,1);
RelErr_pri = zeros(combinations,1);
innit_pri = zeros(combinations,1);
outit_pri = zeros(combinations,1);

time_rspri = zeros(combinations,1);
Res_rspri = zeros(combinations,1);
RelErr_rspri = zeros(combinations,1);
innit_rspri = zeros(combinations,1);
outit_rspri = zeros(combinations,1);

time_pri_tol = zeros(combinations,1);
Res_pri_tol = zeros(combinations,1);
RelErr_pri_tol = zeros(combinations,1);
innit_pri_tol = zeros(combinations,1);
outit_pri_tol = zeros(combinations,1);

time_cgls = zeros(combinations,1);
Res_cgls = zeros(combinations,1);
RelErr_cgls = zeros(combinations,1);
its_cgls = zeros(combinations,1);

time_cgls_projected = zeros(combinations,1);
Res_cgls_projected = zeros(combinations,1);
RelErr_cgls_projected = zeros(combinations,1);
its_cgls_projected = zeros(combinations,1);

time_nnfcgls = zeros(combinations,1);
Res_nnfcgls = zeros(combinations,1);
RelErr_nnfcgls = zeros(combinations,1);
innit_nnfcgls = zeros(combinations,1);
outit_nnfcgls = zeros(combinations,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0; i_tol = 0; cglsId = 0;
for noiseId = 1:length(NoiseLevels)
    for imgId=1:length(true_img)
        clear A; clear b; clear x; clear ProbInfo; clear bn; clear NoiseInfo;
        NoiseLevel = NoiseLevels(noiseId);
        options.trueImage = true_img{imgId};
        [A, b, x, ProbInfo] = PRblurmotion(options);
        [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
        options = IRset('x_true', x, 'RegParam', 'discrep', 'NoStop', 'off', 'NoiseLevel', NoiseLevel);
        
        cglsId = cglsId + 1;
        p_img_cgls{cglsId} = true_img{imgId};
        p_noiseLevel_cgls(cglsId) = NoiseLevel;
        
        % CGLS stopped by DP
        tic;
        [X_cgls_dp, IterInfo_cgls_dp] = IRcgls(A, bn, options);
        time_cgls(cglsId) = toc;
        its_cgls(cglsId) = IterInfo_cgls_dp.its;
        Res_cgls(cglsId) = IterInfo_cgls_dp.Rnrm(end);
        RelErr_cgls(cglsId) = IterInfo_cgls_dp.Enrm(end);

        % MCGLS
        tic;
        [X_cgls_proj, IterInfo_cgls_proj] = IRcgls_projected(A, bn, options);
        time_cgls_projected(cglsId) = toc;
        its_cgls_projected(cglsId) = IterInfo_cgls_proj.its;
        Res_cgls_projected(cglsId) = IterInfo_cgls_proj.Rnrm(end);
        RelErr_cgls_projected(cglsId) = IterInfo_cgls_proj.Enrm(end);

        % NNFCGLS (maxin = 10; maxout = 30 --> like PRI)
        tic;
        [X_nnfcgls, IterInfo_nnfcgls] = IRnnfcgls(A, bn, options);
        time_nnfcgls(cglsId) = toc;
        innit_nnfcgls(cglsId) = IterInfo_nnfcgls.its;
        outit_nnfcgls(cglsId) = IterInfo_nnfcgls.itsInOut(end,1);
        Res_nnfcgls(cglsId) = IterInfo_nnfcgls.Rnrm(end);
        RelErr_nnfcgls(cglsId) = IterInfo_nnfcgls.Enrm(end);
        
        for alphaId = 1:length(alpha)
            for tauId = 1:length(tau)
                i = i+1;
                p_img_pri{i} = true_img{imgId};
                p_noiseLevel_pri(i) = NoiseLevel;
                p_alpha(i) = alpha(alphaId);
                p_tau(i) = tau(tauId);
                
                % PRI
                nx = norm(x);
                [X, rho_R, eta_X, outit_pri(i), innit_pri(i), time_pri(i)] = projected_rest_it(A,bn,p_alpha(i),NoiseLevel,l, maxit, method, p_tau(i));
                RelErr_pri(i) = norm(x-X(:,outit_pri(i)))/nx;
                Res_pri(i) = rho_R(outit_pri(i))/norm(bn);                    
                                    
                % RSPRI
                [X2, rho_R2, eta_X2, outit_rspri(i), innit_rspri(i), time_rspri(i), m] = restrictedPRI(A,bn,p_alpha(i),NoiseLevel,l, maxit,beta, method, p_tau(i));
                RelErr_rspri(i) = norm(x-X2(:,outit_rspri(i)))/nx;
                Res_rspri(i) = rho_R2(outit_rspri(i))/norm(bn);
                
                for tolId = 1:length(tol)
                    i_tol = i_tol+1;
                    p_img_itol{i_tol} = true_img{imgId};
                    p_noiseLevel_itol(i_tol) = NoiseLevel;
                    p_alpha_itol(i_tol) = alpha(alphaId);
                    p_tau_itol(i_tol) = tau(tauId);   
                    p_tol_itol(i_tol) = tol(tolId);
                    
                    % PRI_tol
                    [X_tol, rho_R_tol, eta_X_tol, outit_pri_tol(i_tol), innit_pri_tol(i_tol), time_pri_tol(i_tol)] = pri_tol(A,bn,p_alpha_itol(i_tol),NoiseLevel,l, maxit, method, p_tol_itol(i_tol), p_tau_itol(i_tol));
                    RelErr_pri_tol(i_tol) = norm(x-X_tol(:,outit_pri_tol(i_tol)))/nx;
                    Res_pri_tol(i_tol) = rho_R_tol(outit_pri_tol(i_tol))/norm(bn);

                end
            end
        end
    end
end

p_img_cgls = string(p_img_cgls)';
p_noiseLevel_cgls = p_noiseLevel_cgls';

p_img_pri = string(p_img_pri)';
p_noiseLevel_pri = p_noiseLevel_pri';
p_alpha = p_alpha';
p_tau = p_tau';

p_img_itol = string(p_img_itol)';
p_noiseLevel_itol = p_noiseLevel_itol';
p_alpha_itol = p_alpha_itol';
p_tau_itol = p_tau_itol'; 
p_tol_itol = p_tol_itol'; 

T_cgls = table(p_img_cgls, p_noiseLevel_cgls, RelErr_cgls(1:cglsId), Res_cgls(1:cglsId), its_cgls(1:cglsId), time_cgls(1:cglsId));
writetable(T_cgls,'AnalysisCGLS_motion.csv')

T_cgls_projected = table(p_img_cgls, p_noiseLevel_cgls, RelErr_cgls_projected(1:cglsId), Res_cgls_projected(1:cglsId), its_cgls_projected(1:cglsId), time_cgls_projected(1:cglsId));
writetable(T_cgls_projected,'AnalysisCGLS_projected_motion.csv')

T_nnfcgls = table(p_img_cgls, p_noiseLevel_cgls, RelErr_nnfcgls(1:cglsId), Res_nnfcgls(1:cglsId), outit_nnfcgls(1:cglsId), innit_nnfcgls(1:cglsId), time_nnfcgls(1:cglsId));
writetable(T_nnfcgls,'AnalysisNNFCGLS_motion.csv')

T_pri = table(p_img_pri, p_noiseLevel_pri, p_alpha, p_tau, RelErr_pri(1:i), Res_pri(1:i), outit_pri(1:i), innit_pri(1:i), time_pri(1:i));
writetable(T_pri,'AnalysisPRI_motion.csv')

T_rspri = table(p_img_pri, p_noiseLevel_pri, p_alpha, p_tau, RelErr_rspri(1:i), Res_rspri(1:i), outit_rspri(1:i), innit_rspri(1:i), time_rspri(1:i));
writetable(T_rspri,'AnalysisRSPRI_motion.csv')

T_pri_tol = table(p_img_itol, p_noiseLevel_itol, p_alpha_itol, p_tau_itol, p_tol_itol, RelErr_pri_tol(1:i_tol), Res_pri_tol(1:i_tol), outit_pri_tol(1:i_tol), innit_pri_tol(1:i_tol), time_pri_tol(1:i_tol));
writetable(T_pri_tol,'AnalysisPRI_tol_motion.csv')