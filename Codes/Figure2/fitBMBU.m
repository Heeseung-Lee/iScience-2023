clear all; close all; clc;

savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-9);
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

bnds                = {[10^-5 5],[-5 5],[10^-5 100],[10^-5 5]}; % {sigma_m, mu_0, sigma_0, kappa}
InitRange           = bnds;
stages              = [1 2];
nIter4InitRandomize = 1000;
nBest4SecondStage   = 20;

subjects    = 1;
thrRT       = 0.3;

fitConds.thrRT = thrRT;
for iStage = stages
    fitConds.iStage = iStage;
    switch iStage
        case 1
            stagename                   = 'First';
            fitConds.options            = optimset('MaxFunEvals',50,'MaxIter',50);
            fitConds.nCallfminsearch    = 1;
        case 2
            stagename                   = 'Second';
            fitConds.stagename          = 'SecondStage';
            fitConds.options            = optimset('MaxFunEvals',10^5,'MaxIter',10^5,'TolFun',10^-7,'TolX',10^-7);
            fitConds.nCallfminsearch    = 2;
    end
    switch iStage
        case 1
            fitConds.InitRange           = InitRange;
            nIter                        = nIter4InitRandomize;
            fitConds.nIter               = nIter;
            fitConds.nIter4InitRandomize = nIter;
        case 2
            nIter                       = nBest4SecondStage;
            fitConds.nIter              = nIter;
            fitConds.nBest4SecondStage  = nIter;
    end
    fitConds.bnds           = bnds;
    fitConds.maxTrialLag    = 7;
    idir = ['./' stagename 'Stage/' ];
    jdir = './FirstStage/';
    if isempty(dir(idir)) == 1
        mkdir(idir)
    end
    for iSub = subjects
        imat1 = isempty(dir(['./' idir '/' num2str(iSub) '.mat']));
        imat2 = isempty(dir(['./' idir '/' num2str(iSub) '_computing.mat']));
        if imat1*imat2 == 1
            save([idir '/' num2str(iSub) '_computing.mat'],'iSub')
            if iStage == 2
                load([jdir '/' num2str(iSub) '.mat'],'fitResults')
            end
            fitConds.iSub   = iSub;
            fitConds.Stim   = Stims{iSub};
            fitConds.Chc    = Chcs{iSub};
            fitConds.RT     = RTs{iSub};
            %
            fit_sigma_m     = NaN(nIter,1);
            fit_mu_0        = NaN(nIter,1);
            fit_sigma_0     = NaN(nIter,1);
            fit_kappa       = NaN(nIter,1);
            %
            init_sigma_m    = NaN(nIter,1);
            init_mu_0       = NaN(nIter,1);
            init_sigma_0    = NaN(nIter,1);
            init_kappa      = NaN(nIter,1);
            %
            minus_sum_log_Lh = NaN(nIter,1);
            for iIter = 1:nIter
                fitConds.iIter = iIter;
                %
                if iStage == 2
                    fitConds.BestParams = ...
                        [fitResults.fit_sigma_m(iIter) fitResults.fit_mu_0(iIter) fitResults.fit_sigma_0(iIter)...
                        fitResults.fit_kappa(iIter)];
                end
                [ifitParams,iminus_sum_log_Lh,guessIn] = fit_BMBU(fitConds);
                % parameters =  {sigma_m, mu_0, sigma_0, kappa, LR}
                fit_sigma_m(iIter)      = ifitParams(1);
                fit_mu_0(iIter)         = ifitParams(2);
                fit_sigma_0(iIter)      = ifitParams(3);
                fit_kappa(iIter)        = ifitParams(4);
                %
                init_sigma_m(iIter)     = guessIn(1);
                init_mu_0(iIter)        = guessIn(2);
                init_sigma_0(iIter)     = guessIn(3);
                init_kappa(iIter)       = guessIn(4);

                minus_sum_log_Lh(iIter) = iminus_sum_log_Lh;
            end
            %
            [minus_sum_log_Lh,sInd]     = sort(minus_sum_log_Lh);
            fitResults                  = [];
            fitResults.minus_sum_log_Lh = minus_sum_log_Lh;
            fitResults.fit_sigma_m      = fit_sigma_m(sInd);
            fitResults.fit_mu_0         = fit_mu_0(sInd);
            fitResults.fit_sigma_0      = fit_sigma_0(sInd);
            fitResults.fit_kappa        = fit_kappa(sInd);
            %
            fitConds.init_sigma_m       = init_sigma_m(sInd);
            fitConds.init_mu_0          = init_mu_0(sInd);
            fitConds.init_sigma_0       = init_sigma_0(sInd);
            fitConds.init_kappa         = init_kappa(sInd);
            %
            plot_fittedParams(fitConds,fitResults,idir)
            %
            save([idir '/' num2str(iSub) '.mat'],'fitConds','fitResults')
            delete([idir '/' num2str(iSub) '_computing.mat'])
        end
    end
end