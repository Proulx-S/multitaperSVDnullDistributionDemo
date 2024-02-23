function plotPerm4(svd,f0List)
if ~exist('f0List','var'); f0List = []; end


%% Plot coherence that include 1,2,3...N first modes in the numerator
f = squeeze(svd.fNull);
s = squeeze(svd.s(:,:,ismember(svd.f,f)));
coh     = squeeze(svd.coh(:,:,ismember(svd.f,f)));
cohNull90 = squeeze(svd.cohNull90(:,:,ismember(svd.f,f)));
sNullMean = squeeze(svd.sNullMean);
sNull90 = squeeze(svd.sNull90);
hFig1 = figure('WindowStyle','docked');
hTile = tiledlayout(4,1); hTile.TileSpacing = 'tight'; hTile.Padding = 'tight';
ax1 = nexttile([3 1]);
for i = 1:size(s,1)
    cohNrm(i,:) = sum(coh(1:i,:),1);
    cohNrm(i,:) = cohNrm(i,:)-i/size(s,1);
    cohNrm(i,:) = cohNrm(i,:) * size(s,1) / (size(s,1)-i);

    cohNull90(i,:,:) = sum(cohNull90(1:i,:,:),1);
    cohNull90(i,:,:) = cohNull90(i,:,:)-i/size(s,1);
    cohNull90(i,:,:) = cohNull90(i,:,:) * size(s,1) / (size(s,1)-i);

    plot(f,cohNrm(i,:)); hold on
    label{i} = [num2str(i) '-mode coherence'];
end
grid on
grid minor
axis tight
ylim([0 1])
ylabel('coherence')
if ~isempty(f0List)
    xline(f0List,'g');
end
legend(label)
ax1.XAxis.Visible = 'off';


%% Normalize singular spectra
% s = s./mean(sP(1,:,:),3);
% sP = sP./mean(sP(1,:,:),3);
% s = s./mean(mean(sP,3),1);
% sP = sP./mean(mean(sP,3),1);


%% Compute rank
[~,b] = min(abs(svd.f - f0List),[],3);
sPthresh = sNull90(:,:,2);
% sPthresh = sNullMean;
for i = 1:size(s,2)
    ind = find(s(:,i)<sPthresh(:,i),1,'first');
    if isempty(ind)
        sRank(i) = nan;
    else
        if ind==1
            sRank(i) = 0;
            svAtCrossing = s(1,i);
        else
            sRank(i) = interp1(s(ind-1:ind,i)-sPthresh(ind-1:ind,i),ind-1:ind,0);
            svAtCrossing = interp1(ind-1:ind,s(ind-1:ind,i),sRank(i));
        end
    end

    
    %%% Plot distribution at selected frequencies
    if ismember(i,b)
        f0 = svd.f(i);
        figure('WindowStyle','docked');
        hP1 = plot(squeeze(s(:,i)),'--ok','LineWidth',3);
        % hP1 = plot(squeeze(s(:,i)./s(1,i)),'--ok','LineWidth',3);
        % hP1 = plot(squeeze(s(:,i)./mean(sP(1,i,:),3)),'--ok','LineWidth',3);
        % hP1 = plot(squeeze(s(:,i)./mean(mean(sP(:,i,:),3),1)),'--ok','LineWidth',3);
        hP1.MarkerFaceColor = 'w';
        hold on
        hP2 = plot(squeeze(sNull90(:,i,:)),'r');
        % hP2 = plot(squeeze(sP(:,i,:)./sP(1,i,:)),'r');
        % hP2 = plot(squeeze(sP(:,i,:)./mean(sP(1,i,:),3)),'r');
        % hP2 = plot(squeeze(sP(:,i,:)./mean(mean(sP(:,i,:),3),1)),'r');
        xlabel('rank'); ylabel('singular value'); xticks(1:size(s,1))
        title(num2str(f0,'%0.3fHz'))
        % plot(sPthresh(:,i),'m')
        uistack(hP1,'top')
        if sRank(i)==0
            xline(1)
        else
            xline(sRank(i))
        end
        yline(svAtCrossing)
        legend([hP1(1) hP2(1)],{'original' ['90% of null distribution (' num2str(svd.nPerm) ' random permutations)']})
    end



    %%% Compute rank agnostic coherence
    if isnan(sRank(i))
        num=nan;
        numInterp=nan;
    else
        r = 1:1:size(s,1);
        rInterp = 1:0.001:size(s,1);
        sInterp = interp1(r,s(:,i),rInterp,'linear',s(1,i));
        
        if sRank(i)==0
            num       = s(r<=1,i);
            numInterp = sInterp(rInterp<=1);
        else
            num       = s(r<=floor(sRank(i)),i);
            numInterp = sInterp(rInterp<=sRank(i));
        end
        nom = s(:,i);
        nomInterp = sInterp;
    end
    lowRankCoh(i) = sum(num.^2)./sum(nom.^2);
    lowRankCohInterp(i) = sum(numInterp.^2)./sum(nomInterp.^2);
end



%% Plot rank
figure(hFig1)
ax2 = nexttile([1 1]);
hPlot1 = plot(f,sRank); hold on
hPlot2 = plot(f,floor(sRank),'--');
axis tight
ylim([-0.1 max(sRank)+0.1])
if ~isempty(f0List)
    xline(f0List,'g');
end
% xlim([0 0.35])
grid on
xlabel(hTile,'Hz');
linkaxes([ax1 ax2],'x')
legend([hPlot1 hPlot2],{'crossing point interpolation' 'integer'})
ylabel('rank')
ax2.YTick = ax2.YTick(round(ax2.YTick)==ax2.YTick);
title(hTile,'Coherence with fixed number of modes')






%% Plot rank agnostic coherence
hFig2 = figure('WindowStyle','docked');
hTile = tiledlayout(4,1); hTile.TileSpacing = 'tight'; hTile.Padding = 'tight';
ax1 = nexttile([3 1]);
hPlot1 = plot(f,coh(1,:),'--'); hold on
hPlot2 = plot(f,lowRankCoh(1,:));
hPlot3 = plot(f,lowRankCohInterp(1,:));
% xlim(f([1 end]))
grid on
grid minor
if ~isempty(f0List)
    xline(f0List,'g');
end
ylabel('coherence')
axis tight
ylim([ 0 1])
ax = gca; ax.XAxis.Visible = 'off';
legend([hPlot1 hPlot2 hPlot3],{'1st-mode' 'low-rank' 'low-rank (upsampled)'})
uistack(hPlot1,'top')

ax2 = nexttile([1 1]);
hPlot1 = plot(f,sRank); hold on
hPlot2 = plot(f,floor(sRank),'--');
axis tight
ylim([-0.1 max(sRank)+0.1])
if ~isempty(f0List)
    xline(f0List,'g');
end
% xlim(f([1 end]))
grid on
legend([hPlot1 hPlot2],{'crossing point interpolation' 'integer'})
ylabel('rank')
ax2.YTick = 0:max(ax2.YTick);

xlabel(hTile,'Hz')
linkaxes([ax1 ax2],'x')
title(hTile,'Coherence with adaptive number of modes')

