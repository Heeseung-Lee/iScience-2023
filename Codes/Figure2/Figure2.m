clear all; close all; clc

nSub            = 41;
maxTrialLag     = 5;
savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-9);
cd(savepath)

addpath('../subfunctions/')

%% load Data

fMRIbehav = load('../../Data/raw/fMRI/behavior.mat');
pupilbehav = load('../../Data/raw/pupil/behavior.mat');
load('../../Data/BMBU/FittedParameters.mat');
load('../../Data/BMBU/Choice_LogisticRegression.mat');

Chcs    = [fMRIbehav.choice; pupilbehav.choice];
Stims   = [fMRIbehav.stimulus; pupilbehav.stimulus];
RTs     = [fMRIbehav.rt; pupilbehav.rt];

%% BMBU

pL = cell(nSub,1);
Conds.thrRT          = 0.3;
Conds.maxTrialLag    = 7;
for iSub = 1:nSub
    Conds.iSub       = iSub;
    Conds.Stim       = Stims{iSub};
    Conds.Chc        = Chcs{iSub};
    Conds.RT         = RTs{iSub};
    %
    Conds.sigma_m    = Params.sigma_m(iSub);
    Conds.mu_0       = Params.mu_0(iSub);
    Conds.sigma_0    = Params.sigma_0(iSub);
    Conds.kappa      = Params.kappa(iSub);
    %
    [~,~,pL{iSub}]      = get_LogLik(Conds);
end

%%

figure(1)
clf

ls  = {':','--','-'};
dx  = [-1 0 1]*0;
lw  = 1.2;
ms  = 8;
fz  = 15;
cls = lines(5);
modelcolor = [1 1 1]*0.7;

% proportion of large choice

subplot(1,2,1)
hold on
PL      = NaN(3,3,nSub);
sPL     = NaN(3,3,nSub);
coefpl  = NaN(nSub,4);
for iSub = 1:nSub
    is      = Stims{iSub} + 2;
    id      = Chcs{iSub}/2 + .5;
    irt     = RTs{iSub};
    nT      = size(is,1);
    nR      = size(is,2);
    ps      = [NaN(1,nR); is(1:end-1,:)];
    pd      = [NaN(1,nR); id(1:end-1,:)];
    inan    = isnan(is + id + ps + pd) | (irt < 0.3);
    for i = 1:3
        for p = 1:3
            iInd            = (is == i) & (ps == p) & ~inan;
            PL(i,p,iSub)    = mean(id(iInd),'omitnan');
            sPL(i,p,iSub)   = mean(pL{iSub}(iInd),'omitnan');
        end
    end
    ynan            = isnan(id(:));
    x               = [is(:) ps(:) pd(:)];
    x               = zscore(x(~inan(:),:));
    y               = id(:);
    y               = y(~inan(:));
    coefpl(iSub,:)  = glmfit(x,y,'binomial','link','logit');
end

for i = 1:3

    iPL         = squeeze(sPL(i,:,:));
    m           = mean(iPL,2,'omitnan');
    plot((1:3)+dx(i),m,'LineStyle',ls{i},'color',modelcolor,'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')

    iPL         = squeeze(PL(i,:,:))';
    m           = mean(iPL,'omitnan');
    [~,~,ci]    = ttest(iPL);
    for p = 1:3
        plot([1 1]*p+dx(i),ci(:,p),'k-_','LineWidth',lw)
    end
    plot((1:3)+dx(i),m,'ok','LineStyle',ls{i},'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')

end
ylabel('proportion of large choice')
set(gca,'xtick',1:3,'xticklabel',-1:1,'fontsize',fz,'ytick',0:0.25:1)
ylim([0 1])
xlim([0.7 3.2])
grid on

% logistic regression

subplot(1,2,2)
hold on
coefs       = NaN(nSub,2*maxTrialLag + 1);
for iSub = 1:nSub
    is          = Stims{iSub} + 2;
    id          = Chcs{iSub}/2 + .5;
    irt         = RTs{iSub};
    nT          = size(is,1);
    nR          = size(is,2);
    inan        = isnan(is + id) | (irt < 0.3);
    id(inan)    = NaN;

    js = is(:);
    jc = [];
    for iTB = 1:maxTrialLag
        ps  = [NaN(iTB,nR); is(1:end-iTB,:)];
        pc  = [NaN(iTB,nR); id(1:end-iTB,:)];
        js  = cat(2,js,ps(:));
        jc  = cat(2,jc,pc(:));
    end

    x                   = [js jc];
    iIndFirst           = [ones(1,nR); zeros(nT-1,nR)];
    iIndFirst           = iIndFirst(:)==1;
    x(iIndFirst,:)      = [];
    y                   = id(:);
    y(iIndFirst,:)      = [];

    inan            = sum(isnan(x),2) == size(x,2) | isnan(y);
    x(inan,:)       = NaN;
    x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
    x(isnan(x))     = 0;
    x(inan,:)       = [];
    x               = zscore(x);
    y               = y(~inan);

    jcoef           = glmfit(x,y,"binomial","link","logit");
    coefs(iSub,:)   = jcoef(2:end);
end

m           = mean(coefs);
sm          = mean(sbeta);
[~,p,ci]    = ttest(coefs);
n           = length(m);
cl          = cls(3:5,:);
bar(1,sm(1),'FaceColor','w','EdgeColor',cl(1,:),'LineWidth',lw)
bar(2:(1+maxTrialLag),sm(2:(1+maxTrialLag)),'FaceColor','w','EdgeColor',cl(2,:),'LineWidth',lw)
bar((2+maxTrialLag):n,sm((2+maxTrialLag):n),'FaceColor','w','EdgeColor',cl(3,:),'LineWidth',lw)

plot([1 1],ci(:,1),'_-','Color',cl(1,:),'LineWidth',lw)
for i = 1:maxTrialLag
    plot([1 1]+i,ci(:,1+i),'_-','Color',cl(2,:),'LineWidth',lw)
    plot([1 1]+maxTrialLag+i,ci(:,1+maxTrialLag+i),'_-','Color',cl(3,:),'LineWidth',lw)
end
plot(1,m(1),'o','Color',cl(1,:),'LineWidth',lw,'MarkerFaceColor','w')
plot(2:(1+maxTrialLag),m(2:(1+maxTrialLag)),'o','Color',cl(2,:),'LineWidth',lw,'MarkerFaceColor','w')
plot((2+maxTrialLag):n,m((2+maxTrialLag):n),'o','Color',cl(3,:),'LineWidth',lw,'MarkerFaceColor','w')
idx = 0.32;
idy = 0.15;
for i = 1:n
    ip = p(i);
    if ip<0.05
        if ip<0.001
            ix = [-1 0 1]*idx;
        elseif ip<0.01
            ix = [-1 1]/2*idx;
        else
            ix = 0;
        end
        m = mean(ci(:,i));
        if m > 0
            plot(i+ix,ci(2,i)+idy,'k*','LineWidth',lw,'MarkerSize',ms)
        else
            plot(i+ix,ci(1,i)-idy,'k*','LineWidth',lw,'MarkerSize',ms)
        end
    end
end
set(gca,'xtick',[1 2 maxTrialLag+1 maxTrialLag+2 n],'xticklabel',[0 1 maxTrialLag 1 maxTrialLag],'fontsize',fz,'ytick',-0.5:0.5:1.5)
ylabel('reg. coef.')
ylim([-1.1 2])
xlim([0 n+1])
grid on


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1200 600]/150)
saveas(gcf,'Figure2.png')
