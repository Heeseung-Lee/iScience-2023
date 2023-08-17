function [LogLik,sum_LogLik,pL,pre_Stim,pre_sigma_m,mu_Bnd,sigma_Bnd,uncert_Stim,uncert_Bnd,dof] = get_LogLik(fitConds)

Chc     = fitConds.Chc;
RT      = fitConds.RT;
thrRT   = fitConds.thrRT;
%
[mu_Bnd,sigma_Bnd,uncert_Bnd]   = get_boundary(fitConds);

sigma_m         = fitConds.sigma_m;
mu_0            = fitConds.mu_0;
sigma_0         = fitConds.sigma_0;
Stim            = fitConds.Stim;
Wt2Stim         = sigma_0^2/(sigma_0^2 + sigma_m^2);
Wt2Prior        = sigma_m^2/(sigma_0^2 + sigma_m^2);
pre_Stim        = Stim*Wt2Stim + mu_0*Wt2Prior;
pre_sigma_m     = ones(size(pre_Stim))*sigma_m*Wt2Stim;
uncert_Stim     = sqrt((ones(size(pre_Stim))*sigma_m*Wt2Stim).^2 + sigma_0.^2*Wt2Prior.^2);

pL              = normcdf((pre_Stim-mu_Bnd)./sqrt(pre_sigma_m.^2+sigma_Bnd.^2));
%
LogLik              = NaN(size(Chc));
LogLik(Chc == -1)   = 1-pL(Chc == -1);
LogLik(Chc == 1)    = pL(Chc == 1);
LogLik(RT<thrRT)    = NaN;
LogLik(1,:)         = NaN;
sum_LogLik          = sum(sum(log(LogLik(~isnan(LogLik)))));
dof                 = sum(~isnan(LogLik(:)));

end
