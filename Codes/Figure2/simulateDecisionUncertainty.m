clear all; close all; clc;

savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-30);
cd(savepath)

fMRIbehav = load('../../Data/raw/fMRI/behavior.mat');
pupilbehav = load('../../Data/raw/pupil/behavior.mat');
load('../../Data/BMBU/FittedParameters.mat');
load('../../Data/BMBU/Choice_LogisticRegression.mat');

Chcs    = [fMRIbehav.choice; pupilbehav.choice];
Stims   = [fMRIbehav.stimulus; pupilbehav.stimulus];
RTs     = [fMRIbehav.rt; pupilbehav.rt];
addpath('../subfunctions/')

%%
nIter_Goal      = 10000000;
nIter_Buffer    = nIter_Goal;

thrRT           = 0.3;
fitConds.thrRT  = thrRT;

fitConds.maxTrialLag    = 7;
nSub = 1;

asigma_m            = NaN(nSub,1);
amu_0               = NaN(nSub,1);
asigma_0            = NaN(nSub,1);
akappa              = NaN(nSub,1);
for iSub = 1:nSub
    load(['./SecondStage/' num2str(iSub) '.mat'])
    asigma_m(iSub)      = fitResults.fit_sigma_m(1);
    amu_0(iSub)         = fitResults.fit_mu_0(1);
    asigma_0(iSub)      = fitResults.fit_sigma_0(1);
    akappa(iSub)        = fitResults.fit_kappa(1);
end
%
if isempty(dir('DecisionUncertainty'))
    mkdir('DecisionUncertainty')
end
cd('DecisionUncertainty')

for iSub = 1:nSub
    mat1 = isempty(dir(['Sub' num2str(iSub) '.mat']));
    mat2 = isempty(dir(['Sub' num2str(iSub) '_computing.mat']));
    %         mat1 = 1;mat2 = 1;
    if mat1*mat2 == 1
        save(['Sub' num2str(iSub) '_computing.mat'],'iSub')
        sSub            = num2str(iSub);
        fitConds.iSub   = iSub;
        fitConds.Stim   = Stims{iSub};
        fitConds.Chc    = Chcs{iSub};
        fitConds.RT     = RTs{iSub};
        %
        fitConds.sigma_m    = asigma_m(iSub);
        fitConds.mu_0       = amu_0(iSub);
        fitConds.sigma_0    = asigma_0(iSub);
        fitConds.kappa      = akappa(iSub);
        %
        [~,~,~,uSen,sigma_m,uBnd,sBnd,seSen,seBnd] = get_LogLik(fitConds);
        indNaN          = isnan(Chcs{iSub}) | RTs{iSub} < thrRT | isnan(uBnd);
        uSen(indNaN)    = NaN;
        sigma_m(indNaN) = NaN;
        seSen(indNaN)   = NaN;
        uBnd(indNaN)    = NaN;
        sBnd(indNaN)    = NaN;
        seBnd(indNaN)   = NaN;
        iS              = Stims{iSub};
        iS(indNaN)      = NaN;
        iC              = Chcs{iSub};
        iC(indNaN)      = NaN;
        iC(iC==-1)      = 0;
        %
        nT          = size(uSen,1);
        nR          = size(uSen,2);
        sTR         = num2str((nT-1)*nR);
        simul_U     = NaN(size(uSen,1),size(uSen,2));
        cT      = 1;
        cSample = 1;
        for iT = 2:nT
            for iR = 1:nR
                sT = num2str(cT);
                oD = iC(iT,iR);
                if ~isnan(oD)
                    disp(['iSub=' sSub ', iTrial=' sT '/' sTR])
                    iuSen       = uSen(iT,iR);
                    isigma_m    = sigma_m(iT,iR);
                    iseSen      = seSen(iT,iR);
                    iuBnd       = uBnd(iT,iR);
                    isBnd       = sBnd(iT,iR);
                    iseBnd      = seBnd(iT,iR);

                    mSen = [];
                    mBnd = [];
                    iInd = [];
                    for i = 1:100
                        if length(iInd) < nIter_Goal
                            mSen    = cat(1,mSen,normrnd(iuSen,isigma_m,nIter_Buffer,1));
                            mBnd    = cat(1,mBnd,normrnd(iuBnd,isBnd,nIter_Buffer,1));
                            iInd    = find((mSen > mBnd) == oD);
                            mSen    = mSen(iInd);
                            mBnd    = mBnd(iInd);
                        end
                    end
                    mSen    = mSen(1:nIter_Goal);
                    mBnd    = mBnd(1:nIter_Goal);
                    iU      = normcdf(-abs(mSen-mBnd)./sqrt(iseSen^2+iseBnd^2));

                    simul_U(iT,iR)    = mean(iU);
                end
                cT = cT + 1;
            end
        end
        save(['Sub' num2str(iSub) '.mat'],'simul_U') 
        delete(['Sub' num2str(iSub) '_computing.mat'])
    end
end