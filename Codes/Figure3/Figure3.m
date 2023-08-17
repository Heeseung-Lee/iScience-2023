clear all; close all; clc


fMRIbehav   = load('../../Data/raw/fMRI/behavior.mat');
pupilbehav  = load('../../Data/raw/pupil/behavior.mat');
DUs         = load('../../Data/BMBU/DecisionUncertainty.mat');

Chcs    = [fMRIbehav.choice; pupilbehav.choice];
Stims   = [fMRIbehav.stimulus; pupilbehav.stimulus];
RTs     = [fMRIbehav.rt; pupilbehav.rt];

%% load simulation data

nSub    = 41;
nTB     = 5;


ls  = {':','--','-'};
dx  = [-1 0 1]*0;
lw  = 1.2;
ms  = 8;
fz  = 15;
cls = lines(5);
modelcolor = [1 1 1]*0.7;

figure(1)
clf

DFs = NaN(3,3,nSub);
for iSub = 1:nSub
    iRT         = RTs{iSub}(:);
    iu          = DUs.DU{iSub}(:);
    istm        = Stims{iSub};
    ichc        = Chcs{iSub};
    nR          = size(istm,2);
    nT          = size(istm,1);

    pstm        = [NaN(1,nR); istm(1:end-1,:)];
    icong       = istm(:).*ichc(:);
    pcong       = pstm(:).*ichc(:);
    inan        = isnan(sum([iu iRT icong pcong],2));
    jDF         = zscore(iu(~inan,:));
    icong       = icong(~inan) + 2;
    pcong       = pcong(~inan) + 2;   
    for i = 1:3
        for p = 1:3
            iInd            = (icong == i) & (pcong == p);
            DFs(i,p,iSub)   = mean(jDF(iInd));
        end
    end    
end

cl = lines(5);
cl = cl(3:5,:);

subplot(1,2,1)
hold on
for is = 1:3
    irt         = squeeze(DFs(is,:,:))';
    m           = mean(irt,'omitnan');
    plot((1:3)+dx(is),m,'color','k','LineStyle',ls{is},'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')
end
ylabel('decision uncertainty (z)')
% xlabel('D(t)*S(t-1)')
set(gca,'xtick',1:3,'xticklabel',-1:1,'fontsize',fz,'ytick',-1:1)
ylim([-1.6 2.1])
xlim([0.7 3.2])
grid on

% subplot(6,3,18): regression

coefs       = NaN(nSub,2*nTB + 1);
for iSub = 1:nSub
    iRT         = RTs{iSub}(:);
    iu          = DUs.DU{iSub,1}(:);
    istm        = Stims{iSub};
    ichc        = Chcs{iSub};
    nR          = size(istm,2);
    nT          = size(istm,1);

    is = istm(:).*ichc(:);
    ic = [];
    for iTB = 1:nTB
        ps  = [NaN(iTB,nR); istm(1:end-iTB,:)];
        pc  = [NaN(iTB,nR); ichc(1:end-iTB,:)];
        is  = cat(2,is,ps(:).*ichc(:));
        ic  = cat(2,ic,pc(:).*ichc(:));
    end
    x                   = [is ic];
    iIndFirst           = [ones(1,nR); zeros(nT-1,nR)];
    iIndFirst           = iIndFirst(:)==1;
    x(iIndFirst,:)      = [];
    iu(iIndFirst,:)     = [];
    iRT(iIndFirst,:)    = [];

    inan            = sum(isnan(x),2) == size(x,2) | isnan(iu) | iRT<0.3;
    x(inan,:)       = NaN;
    x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
    x(isnan(x))     = 0;
    x(inan,:)       = [];
    x               = zscore(x);
    y               = zscore(iu(~inan));
    
    jcoef           = glmfit(x,y);
    coefs(iSub,:)   = jcoef(2:end);
end
subplot(1,2,2)
hold on
m = mean(coefs);
n = length(m);
[~,p] = ttest(coefs);
bar(1,m(1),'FaceColor',cl(1,:),'EdgeColor',cl(1,:),'LineWidth',lw)
bar(2:(1+nTB),m(2:(1+nTB)),'FaceColor',cl(2,:),'EdgeColor',cl(2,:),'LineWidth',lw)
bar((2+nTB):n,m((2+nTB):n),'FaceColor',cl(3,:),'EdgeColor',cl(3,:),'LineWidth',lw)
idx = 0.3;
idy = 0.015;
set(gca,'xtick',[1 2 nTB+1 nTB+2 n],'xticklabel',[0 1 nTB 1 nTB],'fontsize',fz,'ytick',-0.5:0.5:0.5)
ylabel('reg. coef.')
% xlabel('    D(t)*S(t-i)                 D(t)*D(t-i)  ')
ylim([-1 0.6])
xlim([0 n+1])
grid on


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1300 600]/150)
saveas(gcf,'Figure3.png')
