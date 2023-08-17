clear all; close all; clc

pupilsize   = load('../../Data/preprocessed/pupil.mat');
pupilbehav  = load('/Volumes/ROOT/CSNL_temp/HL/results/RTPupildACC/Data_and_Codes/Data/raw/pupil/behavior.mat');
DU          = load('../../Data/BMBU/DecisionUncertainty.mat');

%%


nTB             = 5;

timeChc         = [-0.5 2.5];
timeStm         = [-2.2 0];

Cutoff          = [100 0.25];
BufferSec       = 0.2;
dpa_thr         = 50;
ResampleHz      = 10;
EventOnset      = 1;
tEventOnset     = [2.2 2.7];
PreresponseOut  = 0;
TonicOut        = 0;
GazeOut         = 1;

nImgSec         = 500;

frameChc        = timeChc*ResampleHz;
nfChc           = frameChc(2)-frameChc(1) + 1;
frameStm        = timeStm*ResampleHz;
nfStm           = frameStm(2)-frameStm(1) + 1;
nfplot          = nfChc + nfStm;
iframe_son      = 2.2*ResampleHz+1;

missout         = 1;
missperiod      = [2 4];

tXhatChc        = [0.5 0.5]; % from choice onset
iIndXhat        = round(ResampleHz*tXhatChc);
nXhat           = length(iIndXhat);

iHz             = 500/ResampleHz;
nfTotal         = 178200;
nfTrial         = ResampleHz*13.2;
nfx_re          = round(linspace(1,nfTotal,nfTotal/iHz));

timeChc_orig    = [-2 3];
frameChc_orig   = timeChc_orig*ResampleHz;
nfChc_orig      = frameChc_orig(2)-frameChc_orig(1) + 1;

nR              = 6;
nT              = 27;
nSub            = 23;

%%

pa = cell(nSub,1);
for iSub = 1:nSub
    ipa_slock       = pupilsize.signalplots{iSub}{4,2};
    jpa_slock       = ipa_slock(:, round(iframe_son + (frameStm(1):frameStm(2))));

    ipa_clock       = pupilsize.signalplots{iSub}{6,2};
    jpa_clock       = ipa_clock(:, round(-frameChc_orig(1) + (frameChc(1):frameChc(2))));

    ixhat           = mean(jpa_clock(:,round(-frameChc(1) + (iIndXhat(1):iIndXhat(2)) + 1)),2,'omitnan'); % confirm this index rightly indicates the peak time point!!!!
    pa{iSub}        = reshape(ixhat,nT,nR);
end

%%

coefs       = NaN(nSub,2*nTB + 1);
coef_xhat   = NaN(nSub,4);
RTs         = NaN(3,3,nSub);
R2          = NaN(nSub,4);
v_cspspc    = NaN(nSub,6);
v_cspspc_signed = NaN(nSub,6);
areg        = [];
for iSub = 1:nSub
    %     fprintf(['iSub=%d\n'],iSub)
    iRT         = pa{iSub}(:);
    iu          = DU.DU{iSub+18,1}(:);
    istm        = pupilbehav.stimulus{iSub};
    ichc        = pupilbehav.choice{iSub};
    nR          = size(istm,2);

    % Xhat
    pstm        = [NaN(1,nR); istm(1:end-1,:)];
    pchc        = [NaN(1,nR); ichc(1:end-1,:)];
    icong       = istm(:).*ichc(:);
    pcong       = pstm(:).*ichc(:);
    pccong      = pchc(:).*ichc(:);
    inan        = isnan(sum([iu iRT icong pcong pccong],2));
    jRT         = zscore(iRT(~inan,:));
    icong       = icong(~inan) + 2;
    pcong       = pcong(~inan) + 2;
    pccong      = pccong(~inan);
    for i = 1:3
        for p = 1:3
            iInd            = (icong == i) & (pcong == p);
            RTs(i,p,iSub)   = mean(jRT(iInd));
        end
    end
    coef_xhat(iSub,:)   = glmfit(zscore([icong pcong pccong]),jRT);

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
    iRT(iIndFirst,:)    = [];

    inan            = sum(isnan(x),2) == size(x,2) | isnan(iRT);
    x(inan,:)       = NaN;
    x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
    x(isnan(x))     = 0;
    x(inan,:)       = [];
    x               = zscore(x);
    y               = zscore(iRT(~inan));

    jcoef           = glmfit(x,y);
    coefs(iSub,:)   = jcoef(2:end);

    % variance partitioning

    iRT                 = pa{iSub};
    pRT                 = [NaN(1,nR); iRT(1:end-1,:)];
    pstm                = [NaN(1,nR); istm(1:end-1,:)];
    pchc                = [NaN(1,nR); ichc(1:end-1,:)];
    pstm2               = [NaN(2,nR); istm(1:end-2,:)];
    pchc2               = [NaN(2,nR); ichc(1:end-2,:)];
    icong               = istm.*ichc;
    pcong               = pstm.*ichc;
    pccong              = pchc.*ichc;
    pcong2              = pstm2.*ichc;
    pccong2             = pchc2.*ichc;

    nT                  = size(iRT,1);
    linx                = 1:nT;
    linRT               = NaN(size(iRT));
    quaRT               = NaN(size(iRT));
    for iR = 1:nR
        iglm            = glmfit(linx,iRT(:,iR));
        linRT(:,iR)     = glmval(iglm,linx,'identity');
        iglm            = glmfit(linx.^2,iRT(:,iR));
        quaRT(:,iR)     = glmval(iglm,linx.^2,'identity');
    end

    ireg                = [iRT(:) pRT(:) linRT(:) icong(:) pcong(:) pcong2(:) pccong(:) pccong2(:)];
    ivalid              = ~isnan(sum(ireg,2));
    ireg                = zscore(ireg(ivalid,:));
    y                   = ireg(:,1);
    prt                 = ireg(:,2);
    trd                 = ireg(:,3);
    cs                  = ireg(:,4);
    ps                  = ireg(:,5);
    pc                  = ireg(:,7);

    sst                 = sum((y - mean(y)).^2);
    in                  = length(y);

    for icase = 1:6
        switch icase
            case 1
                x           = [prt trd cs ps pc];
            case 2
                x           = [trd cs ps pc];
            case 3
                x           = [prt cs ps pc];
            case 4
                x           = [prt trd ps pc];
            case 5
                x           = [prt trd cs pc];
            case 6
                x           = [prt trd cs ps];
        end
        ifull               = glmfit(x,y);
        ssr                 = sum((y - glmval(ifull,x,"identity")).^2);
        r2_model            = 1 - ssr/sst;
        if icase == 1
            r2_full                 = r2_model;
            v_cspspc(iSub,icase)    = r2_full;
        else
            v_cspspc(iSub,icase)    = r2_full - r2_model;
        end
    end
end
m = mean(coef_xhat);
[~,p,~,t] = ttest(coef_xhat);
se = std(coef_xhat)/sqrt(nSub);

signedR2    = (R2(:,1) - R2(:,2:4)).*((coef_xhat(:,2:4) > 0) - 0.5)*2;
mR2(1:3)    = mean(signedR2);
mR2(4:5)    = mR2(2:3)./mR2(1);

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
sp1 = subplot(1,4,1);
hold on
for is = 1:3
    irt         = squeeze(RTs(is,:,:))';
    m           = mean(irt,'omitnan');
    [~,~,ci]    = ttest(irt);
    for i = 1:3
        plot([1 1]*i+dx(is),ci(:,i),'k-_','LineWidth',lw)
    end
    plot((1:3)+dx(is),m,'ko','LineStyle',ls{is},'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')
end
ylabel('normalized pupil (z)')
% xlabel('D(t)*S(t-1)')
set(gca,'xtick',1:3,'xticklabel',-1:1,'fontsize',fz,'ytick',-0.2:0.2:0.4)
ylim([-0.25 0.62])
xlim([0.7 3.2])
grid on
title('pupil')

sp2 = subplot(1,4,2);
hold on
m = mean(signedR2);
[~,~,ci] = ttest(signedR2);
for i = 1:3
    bar(i,m(i),'FaceColor','w','EdgeColor',cl(i,:),'LineWidth',lw)
    plot([i i],ci(:,i),'_-','LineWidth',lw,'Color',cl(i,:))
end
xlim([0 4])
ylim([-0.05 0.02])
grid on
set(gca,'xtick',1:3,'xticklabel',{'','',''},'fontsize',fz,'ytick',-0.04:0.02:0.02)
ylabel('signed R-squared')

sp3 = subplot(1,4,3);
hold on
m = mean(coefs);
[~,p,ci] = ttest(coefs);
n = length(m);
bar(1,m(1),'FaceColor','w','EdgeColor',cl(1,:),'LineWidth',lw)
bar(2:(1+nTB),m(2:(1+nTB)),'FaceColor','w','EdgeColor',cl(2,:),'LineWidth',lw)
bar((2+nTB):n,m((2+nTB):n),'FaceColor','w','EdgeColor',cl(3,:),'LineWidth',lw)
plot([1 1],ci(:,1),'_-','Color',cl(1,:),'LineWidth',lw)
for i = 1:nTB
    plot([1 1]+i,ci(:,1+i),'_-','Color',cl(2,:),'LineWidth',lw)
    plot([1 1]+nTB+i,ci(:,1+nTB+i),'_-','Color',cl(3,:),'LineWidth',lw)
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
if nTB == 1
    set(gca,'xtick',[1 2 n],'xticklabel',[0 1 nTB],'fontsize',fz,'ytick',-0.2:0.1:0.1)
else
    set(gca,'xtick',[1 2 nTB+1 nTB+2 n],'xticklabel',[0 1 nTB 1 nTB],'fontsize',fz,'ytick',-0.2:0.1:0.1)
end
ylabel('reg. coef.')
% xlabel('    D(t)*S(t-i)                 D(t)*D(t-i)  ')
ylim([-0.3 0.15])
xlim([0 n+1])
grid on

%%

lw  = 1.2;
lws = [5, 5, 5];
cls = {[0.5 0.5 0.5],[0 0 0],[1 0 0]};
ps  = [0.05, 0.01, 0.001];
dx  = 0.25;
fz  = 12;
ms  = 5;
x_fromS     = [-2 -1 0];
x_fromC     = [0 1 2];
xval        = [x_fromS x_fromC];
xs          = [x_fromS*ResampleHz + nfStm(1) nfStm - frameChc(1) + 1 + x_fromC*ResampleHz];
lwsig   = 3;

hlpa        = cell(1,2);
paXhat      = NaN(nSub,2);
coefs       = NaN(nSub,nfplot);
coefs_cong  = NaN(nSub,nfplot,3);
for iSub = 1:nSub
    ipa_slock       = pupilsize.signalplots{iSub}{4,2};
    jpa_slock       = ipa_slock(:, iframe_son + (frameStm(1):frameStm(2)));

    ipa_clock       = pupilsize.signalplots{iSub}{6,2};
    jpa_clock       = ipa_clock(:, -frameChc_orig(1) + (frameChc(1):frameChc(2)));

    ipa             = [jpa_slock jpa_clock];
    ixhat           = mean(jpa_clock(:,-frameChc(1) + (iIndXhat(1):iIndXhat(2)) + 1),2,'omitnan'); % confirm this index rightly indicates the peak time point!!!!
    iu              = DU.DU{18+iSub}(:);
    istm            = pupilbehav.stimulus{iSub};
    ichc            = pupilbehav.choice{iSub};
    nR              = size(istm,2);
    pstm            = [NaN(1,nR); istm(1:end-1,:)];
    pchc            = [NaN(1,nR); ichc(1:end-1,:)];
    icong           = istm(:).*ichc(:);
    pcong           = pstm(:).*ichc(:);
    pccong          = pchc(:).*ichc(:);

    inan            = isnan(sum([iu jpa_slock jpa_clock ixhat icong pcong pccong],2));
    ipa(inan,:)     = [];
    iu(inan,:)      = [];
    ixhat(inan,:)   = [];
    icong           = zscore(icong(~inan));
    pcong           = zscore(pcong(~inan));
    pccong          = zscore(pccong(~inan));


    iInd            = iu<median(iu);
    hlpa{1}(iSub,:) = mean(ipa(iInd,:));
    hlpa{2}(iSub,:) = mean(ipa(~iInd,:));

    paXhat(iSub,1)  = mean(ixhat(iInd));
    paXhat(iSub,2)  = mean(ixhat(~iInd));

    ju  = zscore(iu);
    for jf = 1:nfplot
        jpa             = zscore(ipa(:,jf));
        icoef           = glmfit(ju,jpa);
        coefs(iSub,jf)  = icoef(2);

        icoef                   = glmfit([icong pcong pccong],jpa);
        coefs_cong(iSub,jf,:)   = icoef(2:4);
    end
end

subplot(1,4,4)
hold on
cl = lines(5);
cl = {[cl(3,:);[1 0 0];[0 1 0]],[cl(4,:);[1 0 0];[0 1 0]],[cl(5,:);[1 0 0];[0 1 0]]};
ys = [-0.18 0.14 0.15];
for icong = 1:3
    ipa             = coefs_cong(:,1:nfStm,icong);
    y               = mean(ipa);
    [~,p,ci]        = ttest(ipa);
    e               = (ci(2,:)-ci(1,:))/2;
    lineProps.col   = {cl{icong}(1,:)};
    lineProps.style = '-';
    x               = 1:nfStm;
    mseb(x,y,e,lineProps,1);
    for isig = 1:3
        iInd        = p < ps(isig);
        plot(x(iInd),ys(icong)*ones(1,sum(iInd)),'color',cl{icong}(isig,:),'linewidth',lws(isig))
    end

    ipa             = coefs_cong(:,nfStm + (1:nfChc),icong);
    y               = mean(ipa);
    [~,p,e]         = ttest(ipa);
    e               = (e(2,:)-e(1,:))/2;
    lineProps.col   = {cl{icong}(1,:)};
    lineProps.style = '-';
    x               = nfStm + (1:nfChc);
    mseb(x,y,e,lineProps,1);
    for isig = 1:3
        iInd        = p < ps(isig);
        plot(x(iInd),ys(icong)*ones(1,sum(iInd)),'color',cl{icong}(isig,:),'linewidth',lws(isig))
    end

end
plot([0 nfplot+1],[0 0],'k--')
xlim([0 nfplot+1])
ylim([-0.2 0.15])
ylabel('reg. coef.')
xlabel('time from stimulus')
set(gca,'ytick',-0.1:0.1:0.1,'fontsize',fz,'xtick',xs,'XTickLabel',xval)
grid on
%%

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2500 500]/150)
saveas(gcf,'Figure8.png')