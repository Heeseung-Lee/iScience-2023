function [mu_Bnd,sigma_Bnd,uncert_Bnd] = get_boundary(fitConds)
sigma_m         = fitConds.sigma_m;
STMincre        = fitConds.kappa;
mu_0            = fitConds.mu_0;
sigma_0         = fitConds.sigma_0;

Stim            = fitConds.Stim;
maxTrialLag     = fitConds.maxTrialLag;

%% weight

nT          = size(Stim,1);
nR          = size(Stim,2);
mu_Bnd      = NaN(size(Stim));
sigma_Bnd   = NaN(size(Stim));
uncert_Bnd  = NaN(size(Stim));

%

for iT = 2:nT
    if iT >= 2 && iT <= (maxTrialLag+1)
        pre_sigma_m     = sigma_m*fliplr((1+STMincre).^(1:(iT-1)));
    else
        pre_sigma_m     = sigma_m*fliplr((1+STMincre).^(1:maxTrialLag));
    end
    if iT >= 2 && iT <= (maxTrialLag+1)
        pre_Stim     = Stim(1:iT-1,:);
    else
        pre_Stim     = Stim((iT-maxTrialLag):iT-1,:);
    end
    
    sigmas          = [sigma_0 pre_sigma_m];
    Wts             = sigmas.^-2/sum(sigmas.^-2);
    iuBnd           = Wts*[mu_0*ones(1,nR); pre_Stim];
    isBnd           = sqrt(sum(Wts(2:end).^2.*pre_sigma_m.^2));
    iseBnd          = sqrt(sum(sigmas.^2.*Wts.^2));

    mu_Bnd(iT,:)        = iuBnd;
    sigma_Bnd(iT,:)     = isBnd;
    uncert_Bnd(iT,:)    = iseBnd;
end