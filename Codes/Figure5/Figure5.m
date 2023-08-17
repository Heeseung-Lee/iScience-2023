clear all; clc; close all

savepath = matlab.desktop.editor.getActiveFilename;
savepath = savepath(1:end-9);
cd(savepath)

fMRIbehav   = load('../../Data/raw/fMRI/behavior.mat');
pupilbehav  = load('../../Data/raw/pupil/behavior.mat');
DUs         = load('../../Data/BMBU/DecisionUncertainty.mat');

Chcs    = [fMRIbehav.choice; pupilbehav.choice];
Stims   = [fMRIbehav.stimulus; pupilbehav.stimulus];
RTs     = [fMRIbehav.rt; pupilbehav.rt];
DUs     = DUs.DU;

nSub            = 41;
maxTrialLag     = 5;

%%

eq = cell(1,2);
eq{1} = 'rt0~1+s0+s1*rt1+d1*rt1';
eq{2} = 'rt0~1+s0+s1*rt1+d1*rt1+u0';
varname = {'rt0','s0'};
for iTB = 1:maxTrialLag
    varname{2+iTB}      = ['s' num2str(iTB)];
    varname{2+maxTrialLag+iTB}  = ['d' num2str(iTB)];
    eq{1}  = strcat(eq{1},'+',['s' num2str(iTB)],'+',['d' num2str(iTB)]);
    eq{2}  = strcat(eq{2},'+',['s' num2str(iTB)],'+',['d' num2str(iTB)]);
end
varname{end+1} = 'u0';
varname{end+1} = 'rt1';

acoef   = cell(1,2);
for icond = 1:2
    coefs = [];
    for iSub = 1:nSub
        iRT         = RTs{iSub}(:);
        istm        = Stims{iSub};
        ichc        = Chcs{iSub};
        iu          = DUs{iSub}(:);
        prt         = [NaN(1,size(istm,2)); RTs{iSub}(1:end-1,:)];
        prt         = prt(:);
        pu          = [NaN(1,size(istm,2)); DUs{iSub,1}(1:end-1,:)];
        pu          = pu(:);

        nR          = size(istm,2);
        nT          = size(istm,1);

        is = istm(:).*ichc(:);
        ic = [];
        for iTB = 1:maxTrialLag
            ps  = [NaN(iTB,nR); istm(1:end-iTB,:)];
            pc  = [NaN(iTB,nR); ichc(1:end-iTB,:)];
            is  = cat(2,is,ps(:).*ichc(:));
            ic  = cat(2,ic,pc(:).*ichc(:));
        end

        x                   = [is ic iu];
        iIndFirst           = [ones(1,nR); zeros(nT-1,nR)];
        iIndFirst           = iIndFirst(:)==1;
        x(iIndFirst,:)      = [];
        iu(iIndFirst,:)     = [];
        iRT(iIndFirst,:)    = [];
        prt(iIndFirst,:)    = [];
        pu(iIndFirst,:)     = [];

        inan            = sum(isnan(x),2) == size(x,2) | isnan(iu) | iRT<0.3;
        x(inan,:)       = NaN;
        x               = (x - mean(x,'omitnan'))./std(x,'omitnan');
        x(isnan(x))     = 0;
        x(inan,:)       = [];
        x               = zscore(x);

        y               = zscore(iRT(~inan));
        prt             = (prt - mean(prt,'omitnan'))./std(prt,'omitnan');
        prt(isnan(prt)) = 0;
        prt(inan,:)     = [];
        prt             = zscore(prt);
        ipu             = prt;

        iInd    = ones(size(y))==1;
        ieq     = eq{icond};
        ireg    = [y x ipu];
        ireg    = zscore(ireg(iInd,:));

        reg     = cell2table(num2cell(ireg),'variablenames',varname);
        iglm    = fitglm(reg,ieq);

        coefs(iSub,:)   = iglm.Coefficients.Estimate(2:end);
    end
    acoef{icond} = coefs;
end

%%

figure(1)
clf
cls = lines(5);
lw  = 1.2;
ms  = 8;
fz  = 15;
dx  = [-1 1]*0.1;
st  = {'--','-'};
cl  = lines(5);
cl  = cl(3:5,:);

%

for icond = 1:2
    sp2 = subplot(3,1,icond);
    hold on
    icoef = acoef{1,icond};
    if icond == 1
        icoef = [icoef(:,[1:(1+2*maxTrialLag) ([3 4 2])+2*maxTrialLag])];
    else
        icoef = [icoef(:,[1:(1+2*maxTrialLag) ([4 5 3 2])+2*maxTrialLag])];
    end
    m = mean(icoef);
    [~,p,ci] = ttest(icoef);
    n = length(m);
    bar(1,m(1),'FaceColor','w','EdgeColor',cl(1,:),'LineWidth',lw)
    bar(2:(1+maxTrialLag),m(2:(1+maxTrialLag)),'FaceColor','w','EdgeColor',cl(2,:),'LineWidth',lw)
    bar((2+maxTrialLag):(n-2-icond),m((2+maxTrialLag):(n-2-icond)),'FaceColor','w','EdgeColor',cl(3,:),'LineWidth',lw)
    plot([1 1],ci(:,1),'_-','Color',cl(1,:),'LineWidth',lw)
    for i = 1:maxTrialLag
        plot([1 1]+i,ci(:,1+i),'_-','Color',cl(2,:),'LineWidth',lw)
        plot([1 1]+maxTrialLag+i,ci(:,1+maxTrialLag+i),'_-','Color',cl(3,:),'LineWidth',lw)
    end
    
    if icond == 1
        for i = 1:2
            bar(n-3+i,m(n-3+i),'FaceColor','w','EdgeColor',cl(1+i,:),'LineWidth',lw)
            plot([1 1]*(n-3+i),ci(:,n-3+i),'_-','Color',cl(1+i,:),'LineWidth',lw)
        end
        bar(n,m(n),'FaceColor','w','EdgeColor','k','LineWidth',lw)
        plot([1 1]*(n),ci(:,n),'_-','Color','k','LineWidth',lw)
    else
        for i = 1:2
            bar(n-4+i,m(n-4+i),'FaceColor','w','EdgeColor',cl(1+i,:),'LineWidth',lw)
            plot([1 1]*(n-4+i),ci(:,n-4+i),'_-','Color',cl(1+i,:),'LineWidth',lw)
        end
        
        bar(n-1,m(n-1),'FaceColor','w','EdgeColor','k','LineWidth',lw)
        plot([1 1]*(n-1),ci(:,n-1),'_-','Color','k','LineWidth',lw)
        bar(n,m(n),'FaceColor','w','EdgeColor','r','LineWidth',lw)
        plot([1 1]*(n),ci(:,n),'_-','Color','r','LineWidth',lw)

    end

    idx = 0.3;
    idy = 0.04 + 0.01*(icond-1);
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
    if icond == 1
        set(gca,'xtick',[1 2 maxTrialLag+1 maxTrialLag+2 maxTrialLag+(6:9)],'xticklabel',{'0','1',maxTrialLag','1',maxTrialLag,'S_t_-_1*C_t X RT_t_-_1','C_t_-_1*C_t X RT_t_-_1','RT_t_-_1'},'fontsize',fz,'ytick',-0.2:0.1:0.2)
        ylim([-0.32 0.23])
    else
        set(gca,'xtick',[1 2 maxTrialLag+1 maxTrialLag+2 maxTrialLag+(6:10)],'xticklabel',{'0','1',maxTrialLag','1',maxTrialLag,'S_t_-_1*C_t X RT_t_-_1','C_t_-_1*C_t X RT_t_-_1','RT_t_-_1','u_t'},'fontsize',fz,'ytick',-0.2:0.2:0.5)
        ylim([-0.28 0.5])
    end
    ylabel('reg. coef.')
    xlabel('trial lag (i)')
    xlim([0 15+1])
    grid on
end


%%

regs = cell(nSub,1);
for iSub = 1:nSub

    irt         = RTs{iSub};
    iu          = DUs{iSub};
    istm        = Stims{iSub};
    ichc        = Chcs{iSub};

    nR          = size(istm,2);
    nT          = size(istm,1);

    prt         = [NaN(1,nR); irt(1:end-1,:)];
    pu          = [NaN(1,nR); iu(1:end-1,:)];
    pstm        = [NaN(1,nR); istm(1:end-1,:)];
    pchc        = [NaN(1,nR); ichc(1:end-1,:)];

    iscong      = istm.*ichc;
    pscong      = pstm.*ichc;
    pccong      = pchc.*ichc;

    zrt         = (irt - mean(irt(:),'omitnan'))/std(irt(:),'omitnan');
    zu          = (iu - mean(iu(:),'omitnan'))/std(iu(:),'omitnan');
    
    reg = [zrt(:) zu(:) iscong(:) pscong(:) pccong(:) prt(:)];
    reg0 = reg(~isnan(sum(reg,2)),:);
    reg = zscore(reg0);
    reg = cell2table(num2cell(reg),'VariableNames',{'rt0','u0','s0','s1','c1','rt1'});
    
    iglm = fitglm(reg,'rt0~1+u0');
    iresid = zscore(reg.rt0 - iglm.predict);
    prt = reg.rt1;
    c1 = reg0(:,5);
    ireg = [iresid prt c1];
    regs{iSub} = ireg;
end

%

nBin = 5;
rts = cell(1,2);
x = NaN(nSub,nBin);
for iSub = 1:nSub
    rt0 = regs{iSub}(:,1);
    rt1 = regs{iSub}(:,2);
    c1 = (regs{iSub}(:,3)==1)+1;
    bnd = prctile(rt1,linspace(0,100,nBin+1));    
    for i = 1:2
        for ibin = 1:nBin
            iInd = (c1 == i) & (rt1 >= bnd(ibin) & (rt1 < bnd(ibin+1)));
            rts{i}(iSub,ibin) = mean(rt0(iInd));

            x(iSub,ibin) = mean(rt1(iInd));
        end
    end
end

%%
subplot(3,1,3)
hold on
cl = {'k','k'};
ls = {'--','-','-'};
dx = [-1 1 0]*0.05;
lw = 1.2;
fc = {'w','k'};
mx = mean(x);
ci = cell(1,2);
for i = 1:2
    if i == 3
        irt = mean(cat(3,rts{1},rts{2}),3);
    else
        irt = rts{i};
    end
    m = mean(irt);
    [~,~,ci{i}] = ttest(irt);
    for ibin = 1:nBin
        plot([1 1]*ibin+dx(i),ci{i}(:,ibin),'_-','color',cl{i},'linestyle',ls{i},'linewidth',lw)
    end
    plot((1:nBin)+dx(i),m,'o','color',cl{i},'linestyle',ls{i},'markerfacecolor',fc{i},'linewidth',lw)
end

dy = 0.02;
idx = 0.02;
ms = 8;
idy = 0.1;
for ibin = 1:nBin
    [~,ip] = ttest(rts{1}(:,ibin),rts{2}(:,ibin));

    if ip<0.05
        if ip<0.001
            ix = [-1 0 1]*idx;
        elseif ip<0.01
            ix = [-1 1]/2*idx;
        else
            ix = 0;
        end
        m = max([ci{1}(:,ibin);ci{2}(:,ibin)]);
        plot(ibin+ix,m+idy,'k*','LineWidth',lw,'MarkerSize',ms)
    end
end

set(gca,'xtick',1:5,'ytick',-0.2:0.2:0.4,'fontsize',15)
ylim([-0.35 0.4])
grid on
ylabel('residual RT (z)')
xlabel('previous RT bin')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1200 1800]/150)
saveas(gcf,'Figure5.png')


