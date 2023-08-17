clear all; clc


maxTrialLag     = 5;
savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-9);
cd(savepath)

fMRIbehav   = load('../../Data/raw/fMRI/behavior.mat');
pupilbehav  = load('../../Data/raw/pupil/behavior.mat');

Chcs    = [fMRIbehav.choice; pupilbehav.choice];
Stims   = [fMRIbehav.stimulus; pupilbehav.stimulus];
RTs     = [fMRIbehav.rt; pupilbehav.rt];

nSub = 41;

%%

nIter       = 100;

coefs       = NaN(nSub,2*maxTrialLag + 1);
coefs3var   = NaN(nSub,4);
R2          = NaN(nSub,4);
iRTs        = NaN(3,3,nSub);
v_cspspc    = NaN(nSub,6);
v_cspspc_signed = NaN(nSub,6);
areg        = [];
for iSub = 1:nSub

    iRT         = RTs{iSub}(:);
    istm        = Stims{iSub};
    ichc        = Chcs{iSub};
    nR          = size(istm,2);
    nT          = size(istm,1);

    % RT
    pstm        = [NaN(1,nR); istm(1:end-1,:)];
    pchc        = [NaN(1,nR); ichc(1:end-1,:)];
    prt         = [NaN(1,nR); RTs{iSub}(1:end-1,:)];
    icong       = istm(:).*ichc(:);
    pcong       = pstm(:).*ichc(:);
    pccong      = pchc(:).*ichc(:);
    pRT         = prt(:);
    prt         = prt(:);
    inan        = isnan(sum([iRT icong pcong pccong pRT],2));
    jRT         = zscore(iRT(~inan,:));
    pRT         = zscore(pRT(~inan,:));
    icong       = icong(~inan) + 2;
    pcong       = pcong(~inan) + 2;   
    pccong      = pccong(~inan);
    for i = 1:3
        for p = 1:3
            iInd            = (icong == i) & (pcong == p);
            iRTs(i,p,iSub)   = mean(jRT(iInd));
        end
    end    
    coefs3var(iSub,:)   = glmfit(zscore([icong pcong pccong]),jRT);

    % variance partitioning
    y       = jRT;
    sst     = sum((y - mean(y)).^2);
    for j = 1:4
        switch j 
            case 1
                x = zscore([icong pcong pccong]);
            case 2
                x = zscore([pcong pccong]);
            case 3
                x = zscore([icong pccong]);
            case 4
                x = zscore([icong pcong]);
        end
        icoef       = glmfit(x,y);
        yhat        = glmval(icoef,x,"identity");
        ssr         = sum((y - yhat).^2);
        R2(iSub,j)  = 1 - ssr/sst;
    end
    
    % regression
    is = istm(:).*ichc(:);
    ic = [];
    for iTB = 1:maxTrialLag
        ps  = [NaN(iTB,nR); istm(1:end-iTB,:)];
        pc  = [NaN(iTB,nR); ichc(1:end-iTB,:)];
        is  = cat(2,is,ps(:).*ichc(:));
        ic  = cat(2,ic,pc(:).*ichc(:));
    end
    x                   = [is ic];
    iIndFirst           = [ones(1,nR); zeros(nT-1,nR)];
    iIndFirst           = iIndFirst(:)==1;
    x(iIndFirst,:)      = [];
    iRT(iIndFirst,:)    = [];
    prt(iIndFirst,:)    = [];

    inan            = sum(isnan(x),2) == size(x,2) | isnan(iRT) | iRT<0.3;
    x(inan,:)       = NaN;
    x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
    x(isnan(x))     = 0;
    x(inan,:)       = [];
    x               = zscore(x);
    y               = zscore(iRT(~inan));
    y1              = (prt(~inan) - mean(prt(~inan),'omitnan'))/std(prt(~inan),'omitnan');
    
    jcoef           = glmfit(x,y);
    coefs(iSub,:)   = jcoef(2:end);

end
 
signedR2    = (R2(:,1) - R2(:,2:4)).*((coefs3var(:,2:4) > 0) - 0.5)*2;

%%

ls = {':','--','-'};
dx = [-1 0 1]*0.05;
lw = 1.2;
ms = 8;
fz = 15;

cl = lines(5);
cl = cl(3:5,:);

figure(1)
clf
sp1 = subplot(1,3,1);
hold on
for is = 1:3
    irt         = squeeze(iRTs(is,:,:))';
    m           = mean(irt,'omitnan');
    [~,~,ci]    = ttest(irt);
    for i = 1:3
        plot([1 1]*i+dx(is),ci(:,i),'k-_','LineWidth',lw)
    end
    plot((1:3)+dx(is),m,'ko','LineStyle',ls{is},'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')
end
ylabel('normalized RT (z)')
% xlabel('D(t)*S(t-1)')
set(gca,'xtick',1:3,'xticklabel',-1:1,'fontsize',fz,'ytick',-0.3:0.3:0.6)
ylim([-0.4 0.8])
xlim([0.7 3.2])
title('RT')
grid on

sp2 = subplot(1,3,2);
hold on
m = mean(signedR2);
[~,~,ci] = ttest(signedR2);
for i = 1:3
    bar(i,m(i),'FaceColor','w','EdgeColor',cl(i,:),'LineWidth',lw)
    plot([i i],ci(:,i),'_-','LineWidth',lw,'Color',cl(i,:))
end
xlim([0 4])
ylim([-0.08 0.03])
grid on
set(gca,'xtick',1:3,'xticklabel',{'','',''},'fontsize',fz,'ytick',-0.06:0.02:0.02)
ylabel('signed R-squared')

sp3 = subplot(1,3,3);
hold on
m = mean(coefs);
[~,p,ci] = ttest(coefs);
n = length(m);
bar(1,m(1),'FaceColor','w','EdgeColor',cl(1,:),'LineWidth',lw)
bar(2:(1+maxTrialLag),m(2:(1+maxTrialLag)),'FaceColor','w','EdgeColor',cl(2,:),'LineWidth',lw)
bar((2+maxTrialLag):n,m((2+maxTrialLag):n),'FaceColor','w','EdgeColor',cl(3,:),'LineWidth',lw)
plot([1 1],ci(:,1),'_-','Color',cl(1,:),'LineWidth',lw)
for i = 1:maxTrialLag
    plot([1 1]+i,ci(:,1+i),'_-','Color',cl(2,:),'LineWidth',lw)
    plot([1 1]+maxTrialLag+i,ci(:,1+maxTrialLag+i),'_-','Color',cl(3,:),'LineWidth',lw)
end
idx = 0.3;
idy = 0.025;
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
set(gca,'xtick',[1 2 maxTrialLag+1 maxTrialLag+2 n],'xticklabel',[0 1 maxTrialLag 1 maxTrialLag],'fontsize',fz,'ytick',-0.2:0.1:0.1)
ylabel('reg. coef.')
% xlabel('    D(t)*S(t-i)                 D(t)*D(t-i)  ')
ylim([-0.3 0.18])
xlim([0 n+1])
grid on

sp1.Position = sp1.Position.*[0.6 1.5 0.67 0.8];
sp2.Position = sp2.Position.*[0.87 1.5 0.7 0.8];
sp3.Position = sp3.Position.*[0.85 1.5 1.7 0.8];

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1600 500]/150)
saveas(gcf,'Figure4.png')
