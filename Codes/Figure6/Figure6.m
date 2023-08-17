clear all; close all; clc;

savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-10);
cd(savepath)

ringbehav       = load('../../Data/raw/fine_ring/behavior.mat');
pitchbehav      = load('../../Data/raw/fine_pitch/behavior.mat');
feedbackbehav   = load('../../Data/raw/feedback/behavior.mat');

Stm(:,1)    = ringbehav.stimulus;
Stm(:,2)    = pitchbehav.stimulus;
Chc(:,1)    = ringbehav.choice;
Chc(:,2)    = pitchbehav.choice;
RT(:,1)     = ringbehav.rt;
RT(:,2)     = pitchbehav.rt;
Iden(:,1)   = ringbehav.clas;
Iden(:,2)   = pitchbehav.clas;

%%

featurename     = {'ring','pitch'};
nSub        = 58;
Nk          = 9;
nTB         = 5;
acoefs      = cell(1,2);
acoef_xhat  = cell(1,2);
aR2         = cell(1,2);
aRTs        = cell(1,2);

for ifeature = 1:2

    coefs       = NaN(nSub,2*nTB + 1);
    coef_xhat   = NaN(nSub,4);
    R2          = NaN(nSub,4);
    RTs         = NaN(3,3,nSub);

    for iSub = 1:nSub
        iRT         = RT{iSub,ifeature}(:);
        iIden       = Iden{iSub,ifeature} - (Nk+1)/2;
        istm        = Stm{iSub,ifeature};
        ichc        = Chc{iSub,ifeature};

        nR          = size(istm,2);
        nT          = size(istm,1);

        % Xhat
        pIden       = [NaN(1,nR); iIden(1:end-1,:)];
        pchc        = [NaN(1,nR); ichc(1:end-1,:)];
        icong       = iIden(:).*ichc(:);
        pcong       = pIden(:).*ichc(:);
        pccong      = pchc(:).*ichc(:);
        inan        = isnan(sum([iRT icong pcong pccong],2));
        jRT         = zscore(iRT(~inan,:));
        icong       = icong(~inan) + (Nk+1)/2;
        pcong       = pcong(~inan) + (Nk+1)/2;
        pccong      = pccong(~inan);
        for i = 1:Nk
            for p = 1:Nk
                iInd            = (icong == i) & (pcong == p);
                RTs(i,p,iSub)   = mean(jRT(iInd),'omitnan');
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

        inan            = sum(isnan(x),2) == size(x,2) | isnan(iRT) | iRT<0.3;
        x(inan,:)       = NaN;
        x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
        x(isnan(x))     = 0;
        x(inan,:)       = [];
        x               = zscore(x);
        y               = zscore(iRT(~inan));

        jcoef           = glmfit(x,y);
        coefs(iSub,:)   = jcoef(2:end);

    end

    acoefs{ifeature}      = coefs;
    acoef_xhat{ifeature}  = coef_xhat;
    aR2{ifeature}         = R2;
    aRTs{ifeature}        = RTs;

end

%%

figure(1)
clf

mr2 = NaN(2,5);
mss = NaN(2,4);
ps = NaN(2,4);
for ifeature = 1:2

    signedR2        = (aR2{ifeature}(:,1) - aR2{ifeature}(:,2:4)).*((acoef_xhat{ifeature}(:,2:4) > 0) - 0.5)*2;
    imr2            = mean(signedR2);
    imr2(4:5)       = imr2(2:3)/imr2(1);
    mr2(ifeature,:) = imr2;

    ls = {':','--','-'};
    dx = [-1 0 1]*0.05;
    lw = 1.2;
    ms = 5;
    fz = 15;

    cl = lines(5);
    cl = cl(3:5,:);
    cls = jet(Nk-2);

    sp1 = subplot(3,3,3*(ifeature-1)+1);
    hold on
    c = 1;
    for is = 3:Nk
        irt         = squeeze(aRTs{ifeature}(is,:,:))';
        m           = mean(irt,'omitnan');
        [~,~,ci]    = ttest(irt);
        plot(1:Nk,m,'o-','color',cls(c,:),'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')
        c = c + 1;
    end
    mss(ifeature,:) = mean(acoef_xhat{ifeature});
    [~,ps(ifeature,:)] = ttest(acoef_xhat{ifeature});
    ylabel('normalized RT (z)')
    set(gca,'xtick',1:2:Nk,'xticklabel',(1:2:Nk)-(Nk+1)/2,'fontsize',fz)
    title(featurename{ifeature})
    grid on

    sp2 = subplot(3,3,3*(ifeature-1)+2);
    hold on
    m = mean(signedR2);
    [~,~,ci] = ttest(signedR2);
    for i = 1:3
        bar(i,m(i),'FaceColor','w','EdgeColor',cl(i,:),'LineWidth',lw)
        plot([i i],ci(:,i),'_-','LineWidth',lw,'Color',cl(i,:))
    end
    xlim([0 4])
    grid on
    set(gca,'xtick',1:3,'xticklabel',{'','',''},'fontsize',fz)
    ylabel('signed R-squared')

    sp3 = subplot(3,3,3*(ifeature-1)+3);
    hold on
    m = mean(acoefs{ifeature});
    [~,p,ci] = ttest(acoefs{ifeature});
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
    set(gca,'xtick',[1 2 nTB+1 nTB+2 n],'xticklabel',[0 1 nTB 1 nTB],'fontsize',fz)
    ylabel('reg. coef.')
    xlim([0 n+1])
    grid on
end


%%
nSub = 30;
coefs       = NaN(nSub,2*nTB + 1);
coef_xhat   = NaN(nSub,4);
R2          = NaN(nSub,4);
RTs         = NaN(5,5,nSub);
v_cspspc    = NaN(nSub,6);
v_cspspc_signed = NaN(nSub,6);
areg        = [];
for iSub = 1:nSub
    iRT         = feedbackbehav.rt{iSub}(:);
    istm        = feedbackbehav.stimulus{iSub};
    ichc        = feedbackbehav.choice{iSub};
    nR          = size(istm,2);
    nT          = size(istm,1);

    % Xhat
    pstm        = [NaN(1,nR); istm(1:end-1,:)];
    pchc        = [NaN(1,nR); ichc(1:end-1,:)];
    icong       = istm(:).*ichc(:);
    pcong       = pstm(:).*ichc(:);
    pccong      = pchc(:).*ichc(:);
    inan        = isnan(sum([iRT icong pcong pccong],2));
    jRT         = zscore(iRT(~inan,:));
    icong       = icong(~inan) + 3;
    pcong       = pcong(~inan) + 3;   
    pccong      = pccong(~inan);
    for i = 1:5
        for p = 1:5
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

    inan            = sum(isnan(x),2) == size(x,2) | isnan(iRT) | iRT<0.3;
    x(inan,:)       = NaN;
    x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
    x(isnan(x))     = 0;
    x(inan,:)       = [];
    x               = zscore(x);
    y               = zscore(iRT(~inan));
    
    jcoef           = glmfit(x,y);
    coefs(iSub,:)   = jcoef(2:end);

end
 
signedR2    = (R2(:,1) - R2(:,2:4)).*((coef_xhat(:,2:4) > 0) - 0.5)*2;

%%

cls = jet(5);
lc = {'k',[1 1 1]*0.6,'k',[1 1 1]*0.6,'k'};
dx = [-2 -1 0 1 2]*0.1;
lw = 1.2;
ms = 5;
fz = 15;

cl = lines(5);
cl = cl(3:5,:);

subplot(3,3,7);
hold on
for is = 1:5
    irt         = squeeze(RTs(is,:,:))';
    m           = mean(irt,'omitnan');
%     [~,~,ci]    = ttest(irt);
%     for i = 1:5
%         plot([1 1]*i+dx(is),ci(:,i),'-_','color',cls(is,:),'LineWidth',lw)
%     end
    plot((1:5)+dx(is),m,'-o','color',cls(is,:),'LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','w')
end
ylabel('normalized RT (z)')
% xlabel('D(t)*S(t-1)')
set(gca,'xtick',1:5,'xticklabel',-2:2,'fontsize',fz,'ytick',-0.2:0.2:0.2)
% ylim([-0.3 0.45])
xlim([0.5 5.5])
title('feedback')
grid on

subplot(3,3,8);
hold on
m = mean(signedR2);
[~,~,ci] = ttest(signedR2);
for i = 1:3
    bar(i,m(i),'FaceColor','w','EdgeColor',cl(i,:),'LineWidth',lw)
    plot([i i],ci(:,i),'_-','LineWidth',lw,'Color',cl(i,:))
end
xlim([0 4])
ylim([-0.04 0.01])
grid on
set(gca,'xtick',1:3,'xticklabel',{'','',''},'fontsize',fz,'ytick',-0.06:0.02:0.02)
ylabel('signed R-squared')

subplot(3,3,9);
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
set(gca,'xtick',[1 2 nTB+1 nTB+2 n],'xticklabel',[0 1 nTB 1 nTB],'fontsize',fz,'ytick',-0.2:0.1:0.1)
ylabel('reg. coef.')
% xlabel('    D(t)*S(t-i)                 D(t)*D(t-i)  ')
ylim([-0.21 0.13])
xlim([0 n+1])
grid on

%%

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2000 1500]/150)
saveas(gcf,'Figure6.png')