% calculate speed vs neighbr distance, directional correlation, and
% radial distribution functions

% issues/to-do:
% - seperate into individual functions for each statistic?

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

% define functions for grpstats
mad1 = @(x) mad(x,1); % median absolute deviation
% alternatively could use boxplot-style confidence intervals on the mean,
% which are 1.57*iqr/sqrt(n) - unclear how justified this is
iqrci = @(x) 1.57*iqr(x)/sqrt(numel(x));
% or one could use a bootstrapped confidence interval
bootserr = @(x) bootci(1e1,{@nanmedian,x},'alpha',0.05,'Options',struct('UseParallel',true));


%% set parameters
dataset = 2;  % enter 1 or 2 to specify which dataset to run the script for
phase = 'stationary'; % 'fullMovie' or 'stationary'
plotDiagnostics = false; % true or false

if dataset ==1
    strains = {'npr1'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'npr1'}; %{'npr1','N2'}
end
wormnums = {'40'};%{'40','HD'};
nStrains = length(strains);
plotColors = lines(nStrains);

if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize = 1e4;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
if plotDiagnostics, visitfreqFig = figure; hold on, end

distBinWidth = 35; % in units of micrometers
maxDist = 4000;
distBins = 0:distBinWidth:maxDist;
dircorrxticks = 0:500:2000;

%% go through strains, densities, movies
for wormnum = wormnums
    speedFig = figure; hold on
    dircorrFig = figure; hold on
    velcorrFig = figure; hold on
    poscorrFig = figure; hold on
    lineHandles = NaN(nStrains,1);
    for strainCtr = 1:nStrains
        %% load data
        if dataset == 1
            [lastFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum{1} '_list.xlsx'],1,'A1:B15','basic');
        elseif dataset == 2
            [lastFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum{1} '_g_list.xlsx'],1,'A1:B15','basic');
        end
        numFiles = length(filenames);
        if wormnum{1} == '40'
            visitfreq = cell(numFiles,1);
        end
        speeds = cell(numFiles,1);
        dxcorr = cell(numFiles,1); % for calculating directional cross-correlation
        vxcorr = cell(numFiles,1); % for calculating velocity cross-correlation
        pairdist = cell(numFiles,1);
        nNbrDist= cell(numFiles,1);
        gr =cell(numFiles,1);
        parfor fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            assert(size(skelData,1)==2)
            assert(size(skelData,2)==2)
            assert(size(skelData,3)==length(trajData.frame_number));
            if all(isnan(skelData(:)))
                warning(['all skeleton are NaN for ' filename])
            end
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));

            if strcmp(phase,'fullMovie')
                lastFrame = numel(unique(trajData.frame_number));
                numFrames = round(lastFrame/frameRate/3);
                framesAnalyzed = randperm(lastFrame,numFrames); % randomly sample frames without replacement
            elseif strcmp(phase,'stationary')
                lastFrame = lastFrames(fileCtr);
                firstFrame = double(round(max(trajData.frame_number)/10)); % cut out the first 10 percent of the movie for stationary phase restriction
                numFrames = round((lastFrame-firstFrame)/frameRate/3);
                framesAnalyzed = randperm((lastFrame-firstFrame),numFrames) + firstFrame; % randomly sample frames without replacement
            end

            %% filter worms
            if plotDiagnostics
                plotIntensitySizeFilter(blobFeats,pixelsize,...
                    intensityThresholds(wormnum{1}),maxBlobSize,...
                    [wormnum{1} ' ' strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
            end
            trajData.has_skeleton = squeeze(~any(any(isnan(skelData)))); % reset skeleton flag for pharynx data
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum{1}),maxBlobSize)...
                &trajData.has_skeleton;

            if strcmp(phase,'stationary')
                phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
                trajData.filtered(~phaseFrameLogInd)=false;
            end

            %% calculate stats
            if wormnum{1} == '40'
                OverallArea = pi*(8300/2)^2;
            else
                OverallArea = peak2peak(trajData.coord_x(trajData.filtered)).*...
                    peak2peak(trajData.coord_y(trajData.filtered)).*pixelsize.^2
            end
            speeds{fileCtr} = cell(numFrames,1);
            dxcorr{fileCtr} = cell(numFrames,1); % for calculating directional cross-correlation
            vxcorr{fileCtr} = cell(numFrames,1); % for calculating velocity cross-correlation
            pairdist{fileCtr} = cell(numFrames,1);
            nNbrDist{fileCtr}= cell(numFrames,1);
            gr{fileCtr} = NaN(length(distBins) - 1,numFrames);
            for frameCtr = 1:numFrames % one may be able to vectorise this
                frame = framesAnalyzed(frameCtr);
                [x ,y] = getWormPositions(trajData, frame, true);
                N = length(x);
                if N>1 % need at least two worms in frame
                    frameLogInd = trajData.frame_number==frame&trajData.filtered;
                    vx = double(blobFeats.velocity_x(frameLogInd));
                    vy = double(blobFeats.velocity_y(frameLogInd));
                    speeds{fileCtr}{frameCtr} = sqrt(vx.^2 + vy.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
                    ox = double(squeeze(skelData(1,1,frameLogInd) - skelData(1,2,frameLogInd)));
                    oy = double(squeeze(skelData(2,1,frameLogInd) - skelData(2,2,frameLogInd)));
                    dxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(ox,oy,true,true); % directional correlation
                    vxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(vx,vy,true,false); % velocity correlation
                    pairdist{fileCtr}{frameCtr} = pdist([x y]).*pixelsize; % distance between all pairs, in micrometer
                    gr{fileCtr}(:,frameCtr) = histcounts(pairdist{fileCtr}{frameCtr},distBins,'Normalization','count'); % radial distribution function
                    gr{fileCtr}(:,frameCtr) = gr{fileCtr}(:,frameCtr)'.*OverallArea ...
                        ./(2*pi*distBins(2:end)*distBinWidth*N*(N-1)/2); % normalisation by N(N-1)/2 as pdist doesn't double-count pairs
                    D = squareform(pairdist{fileCtr}{frameCtr}); % distance of every worm to every other
                    nNbrDist{fileCtr}{frameCtr} = min(D + max(max(D))*eye(size(D)));
                    if (numel(speeds{fileCtr}{frameCtr})~=numel(nNbrDist{fileCtr}{frameCtr}))||(numel(dxcorr{fileCtr}{frameCtr})~=numel(pairdist{fileCtr}{frameCtr}))
                        error(['Inconsistent number of variables in frame ' num2str(frame) ' of ' filename ])
                    end
                end
            end
            % pool data from frames
            speeds{fileCtr} = vertcat(speeds{fileCtr}{:});
            dxcorr{fileCtr} = horzcat(dxcorr{fileCtr}{:});
            vxcorr{fileCtr} = horzcat(vxcorr{fileCtr}{:});
            pairdist{fileCtr} = horzcat(pairdist{fileCtr}{:});
            nNbrDist{fileCtr} = horzcat(nNbrDist{fileCtr}{:});
            % heat map of sites visited - this only makes sense for 40 worm
            % dataset where we don't move the camera
            if strcmp(wormnum{1},'40')&& plotDiagnostics
                siteVisitFig = figure;
                h=histogram2(trajData.coord_x*pixelsize/1000,trajData.coord_y*pixelsize/1000,...
                    'DisplayStyle','tile','EdgeColor','none','Normalization','pdf');
                visitfreq{fileCtr} = h.Values(:);
                cb = colorbar; cb.Label.String = '# visited';
                axis equal
                xlabel('x (mm)'), ylabel('y (mm)')
                title([strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
                set(siteVisitFig,'PaperUnits','centimeters')
                figurename = ['figures/individualRecordings/' strains{strainCtr}...
                    '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited' '_' phase '_data' num2str(dataset)];
                exportfig(siteVisitFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
        %% combine data from multiple files
        nNbrDist = horzcat(nNbrDist{:});
        speeds = vertcat(speeds{:});
        pairdist = horzcat(pairdist{:});
        dxcorr = horzcat(dxcorr{:});
        vxcorr = horzcat(vxcorr{:});
        % bin distance data
        [nNbrDistcounts,nNbrDistBins,nNbrDistbinIdx]  = histcounts(nNbrDist,...
            'BinWidth',distBinWidth,'BinLimits',[min(nNbrDist) maxDist]);
        [pairDistcounts,pairDistBins,pairDistbinIdx]  = histcounts(pairdist,...
            'BinWidth',distBinWidth,'BinLimits',[min(pairdist) maxDist]);
        % convert bin edges to centres (for plotting)
        nNbrDistBins = double(nNbrDistBins(1:end-1) + diff(nNbrDistBins)/2);
        pairDistBins = double(pairDistBins(1:end-1) + diff(pairDistBins)/2);
        % ignore larger distance values and bins with only one element, as this will cause bootsci to fault
        nNdistkeepIdcs = nNbrDistbinIdx>0&ismember(nNbrDistbinIdx,find(nNbrDistcounts>1));
        nNbrDistBins = nNbrDistBins(nNbrDistcounts>1);
        speeds = speeds(nNdistkeepIdcs);
        nNbrDistbinIdx = nNbrDistbinIdx(nNdistkeepIdcs);
        pdistkeepIdcs = pairDistbinIdx>0&ismember(pairDistbinIdx,find(pairDistcounts>1));
        pairDistBins = pairDistBins(pairDistcounts>1);
        dxcorr = dxcorr(pdistkeepIdcs);
        vxcorr = vxcorr(pdistkeepIdcs);
        pairDistbinIdx = pairDistbinIdx(pdistkeepIdcs);
        
        [s_med,s_ci] = grpstats(speeds,nNbrDistbinIdx,{@median,bootserr});
        [corr_o_med,corr_o_ci] = grpstats(dxcorr,pairDistbinIdx,{@median,bootserr});
        [corr_v_med,corr_v_ci] = grpstats(vxcorr,pairDistbinIdx,{@median,bootserr});
        %% plot data
        [lineHandles(strainCtr), ~] = boundedline(nNbrDistBins,smooth(s_med),...
            [smooth(s_med - s_ci(:,1)), smooth(s_ci(:,2) - s_med)],...
            'alpha',speedFig.Children,'cmap',plotColors(strainCtr,:));
        boundedline(pairDistBins,smooth(corr_o_med),[smooth(corr_o_med - corr_o_ci(:,1)), smooth(corr_o_ci(:,2) - corr_o_med)],...
            'alpha',dircorrFig.Children,'cmap',plotColors(strainCtr,:))
        boundedline(pairDistBins,smooth(corr_v_med),[smooth(corr_v_med - corr_v_ci(:,1)), smooth(corr_v_ci(:,2) - corr_v_med)],...
            'alpha',velcorrFig.Children,'cmap',plotColors(strainCtr,:))
        gr = cat(2,gr{:});
        boundedline(distBins(2:end)-distBinWidth/2,nanmean(gr,2),...
            [nanstd(gr,0,2) nanstd(gr,0,2)]./sqrt(nnz(sum(~isnan(gr),2))),...
            'alpha',poscorrFig.Children,'cmap',plotColors(strainCtr,:))
        if  strcmp(wormnum{1},'40')&& plotDiagnostics
            histogram(visitfreqFig.Children,vertcat(visitfreq{:}),'DisplayStyle','stairs','Normalization','pdf')
        end
    end
    %% format and export figures
    for figHandle = [speedFig, dircorrFig, velcorrFig poscorrFig] % common formating for both figures
        set(figHandle,'PaperUnits','centimeters')
    end
    %
    speedFig.Children.YLim = [0 400];
    speedFig.Children.XLim = [0 2000];
    speedFig.Children.XTick = 0:500:2000;
    speedFig.Children.Box = 'on';
    speedFig.Children.XDir = 'reverse';
    ylabel(speedFig.Children,'speed (μm/s)')
    xlabel(speedFig.Children,'distance to nearest neighbour (μm)')
    legend(speedFig.Children,lineHandles,strains)
    figurename = ['figures/correlation/phaseSpecific/speedvsneighbrdistance_' wormnum{1} '_' phase '_data' num2str(dataset)];
    exportfig(speedFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    %
    %         dircorrFig.Children.YLim = [-1 1];
% %     dircorrFig.Children.XLim = [0 2000];
    set(dircorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(dircorrFig.Children,'orientational correlation')
    xlabel(dircorrFig.Children,'distance between pair (μm)')
    legend(dircorrFig.Children,lineHandles,strains)
    figurename = ['figures/correlation/phaseSpecific/dircrosscorr_' wormnum{1} '_' phase '_data' num2str(dataset)];
    exportfig(dircorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    %
    %         velcorrFig.Children.YLim = [-1 1];
% %     velcorrFig.Children.XLim = [0 2000];
    set(velcorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(velcorrFig.Children,'velocity correlation')
    xlabel(velcorrFig.Children,'distance between pair (μm)')
    legend(velcorrFig.Children,lineHandles,strains)
    figurename = ['figures/correlation/phaseSpecific/velcrosscorr_' wormnum{1} '_' phase '_data' num2str(dataset)];
    exportfig(velcorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    %
    poscorrFig.Children.YLim(1) = 0;
    poscorrFig.Children.XLim = [0 2000];
    poscorrFig.Children.XTick = 0:500:2000;
    poscorrFig.Children.YTick = 0:round(poscorrFig.Children.YLim(2));
    poscorrFig.Children.Box = 'on';
    ylabel(poscorrFig.Children,'positional correlation g(r)')
    xlabel(poscorrFig.Children,'distance r (μm)')
    legend(poscorrFig.Children,lineHandles,strains)
    figurename = ['figures/correlation/phaseSpecific/radialdistributionfunction_' wormnum{1} '_' phase '_data' num2str(dataset)];
    exportfig(poscorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    
    if  strcmp(wormnum{1},'40')&& plotDiagnostics
        visitfreqFig.Children.XScale = 'log';
        visitfreqFig.Children.YScale = 'log';
        %         visitfreqFig.Children.XLim = [4e-5 1e-1];
        xlabel(visitfreqFig.Children,'site visit frequency, f')
        ylabel(visitfreqFig.Children,'pdf p(f)')
        legend(visitfreqFig.Children,strains)
        figurename = ['figures/visitfreq_' wormnum{1}];
        exportfig(visitfreqFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end
