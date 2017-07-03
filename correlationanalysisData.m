% calculate speed vs neighbr distance, directional correlation, and
% radial distribution functions

% issues/to-do:
% - seperate into individual functions for each statistic?
% - calculate correlation functions separately for worms in/out of clusters
% (the only non-obvious result may be a difference between clustered worms
% in social vs asocial strains)
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
bootserr = @(x) bootci(1e2,{@nanmedian,x},'alpha',0.05,'Options',struct('UseParallel',true));

distBinwidth = 25; % in units of micrometers
maxDist = 2000;
minDist = 50;
distBins = 0:distBinwidth:maxDist;
dircorrxticks = 0:500:2000;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'npr1','N2'};
nStrains = length(strains);
plotColors = lines(nStrains);
wormnums = {'HD'};%{'40','HD'};
intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize = 1e4;
plotDiagnostics = false;
if plotDiagnostics, visitfreqFig = figure; hold on, end
for wormnum = wormnums
    speedFig = figure; hold on
    dircorrFig = figure; hold on
    velcorrFig = figure; hold on
    poscorrFig = figure; hold on
    lineHandles = NaN(nStrains,1);
    for strainCtr = 1:nStrains
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_list.txt']);
        numFiles = length(filenames);
        if wormnum{1} == '40'
            visitfreq = cell(numFiles,1);
        end
        speeds = cell(numFiles,1);
        dxcorr = cell(numFiles,1); % for calculating directional cross-correlation
        vxcorr = cell(numFiles,1); % for calculating velocity cross-correlation
        pairdist = cell(numFiles,1);
        mindist= cell(numFiles,1);
        gr =cell(numFiles,1);
        for fileCtr = 1:numFiles
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
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate/3);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
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
            %% calculate stats
            speeds{fileCtr} = cell(numFrames,1);
            dxcorr{fileCtr} = cell(numFrames,1); % for calculating directional cross-correlation
            vxcorr{fileCtr} = cell(numFrames,1); % for calculating velocity cross-correlation
            pairdist{fileCtr} = cell(numFrames,1);
            mindist{fileCtr}= cell(numFrames,1);
            gr{fileCtr} = NaN(length(distBins) - 1,numFrames);
            for frameCtr = 1:numFrames % one may be able to vectorise this
                frame = framesAnalyzed(frameCtr);
                [x ,y] = getWormPositions(trajData, frame, true);
                if numel(x)>1 % need at least two worms in frame
                    frameLogInd = trajData.frame_number==frame&trajData.filtered;
                    vx = double(blobFeats.velocity_x(frameLogInd));
                    vy = double(blobFeats.velocity_y(frameLogInd));
                    speeds{fileCtr}{frameCtr} = sqrt(vx.^2 + vy.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
                    ox = double(squeeze(skelData(1,1,frameLogInd) - skelData(1,2,frameLogInd)));
                    oy = double(squeeze(skelData(2,1,frameLogInd) - skelData(2,2,frameLogInd)));
                    dxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(ox,oy,true); % directional correlation
                    vxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(vx,vy,true); % velocity correlation
                    pairdist{fileCtr}{frameCtr} = pdist([x y]).*pixelsize; % distance between all pairs, in micrometer
                    gr{fileCtr}(:,frameCtr) = histcounts(pairdist{fileCtr}{frameCtr},distBins,'Normalization','count'); % radial distribution function
                    gr{fileCtr}(:,frameCtr) = gr{fileCtr}(:,frameCtr)'.*maxDist^2./(2*distBins(2:end)*distBinwidth)...
                        ./numel(pairdist{fileCtr}{frameCtr})*2; % normalisation by N(N-1)
                    D = squareform(pairdist{fileCtr}{frameCtr}); % distance of every worm to every other
                    mindist{fileCtr}{frameCtr} = min(D + max(max(D))*eye(size(D)));
                    if (numel(speeds{fileCtr}{frameCtr})~=numel(mindist{fileCtr}{frameCtr}))||(numel(dxcorr{fileCtr}{frameCtr})~=numel(pairdist{fileCtr}{frameCtr}))
                        error(['Inconsistent number of variables in frame ' num2str(frame) ' of ' filename ])
                    end
                end
            end
            % pool data from frames
            speeds{fileCtr} = vertcat(speeds{fileCtr}{:});
            dxcorr{fileCtr} = horzcat(dxcorr{fileCtr}{:});
            vxcorr{fileCtr} = horzcat(vxcorr{fileCtr}{:});
            pairdist{fileCtr} = horzcat(pairdist{fileCtr}{:});
            mindist{fileCtr} = horzcat(mindist{fileCtr}{:});
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
                    '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited'];
                exportfig(siteVisitFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
        %% combine data from multiple files
        mindistVals = quant(horzcat(mindist{:}),distBinwidth);
        speeds = vertcat(speeds{:});
        pairdistVals = quant(horzcat(pairdist{:}),distBinwidth);
        dxcorr = horzcat(dxcorr{:});
        vxcorr = horzcat(vxcorr{:});
        % ignore long distance data (and short distance which can give
        % error in the calculation for too few samples
        speeds = speeds(mindistVals<=maxDist);
        mindistVals = mindistVals(mindistVals<=maxDist);
        pDistKeepIdcs = pairdistVals<=maxDist&pairdistVals>=minDist;
        dxcorr = dxcorr(pDistKeepIdcs);
        vxcorr = vxcorr(pDistKeepIdcs);
        pairdistVals = pairdistVals(pDistKeepIdcs);
        % bootstrapping will yield an error if any bin has n=1
        [s_med,s_ci] = grpstats(speeds,mindistVals,{@median,bootserr});
        [corr_o_med,corr_o_ci] = grpstats(dxcorr,pairdistVals,{@median,bootserr});
        [corr_v_med,corr_v_ci] = grpstats(vxcorr,pairdistVals,{@median,bootserr});
        %% plot data
        mindistBins = double(unique(mindistVals));
        [lineHandles(strainCtr), ~] = boundedline(mindistBins,smooth(s_med),...
            [smooth(s_med - s_ci(:,1)), smooth(s_ci(:,2) - s_med)],...
            'alpha',speedFig.Children,'cmap',plotColors(strainCtr,:));
        pairdistBins = double(unique(pairdistVals));
        boundedline(pairdistBins,smooth(corr_o_med),[smooth(corr_o_med - corr_o_ci(:,1)), smooth(corr_o_ci(:,2) - corr_o_med)],...
            'alpha',dircorrFig.Children,'cmap',plotColors(strainCtr,:))
        boundedline(pairdistBins,smooth(corr_v_med),[smooth(corr_v_med - corr_v_ci(:,1)), smooth(corr_v_ci(:,2) - corr_v_med)],...
            'alpha',velcorrFig.Children,'cmap',plotColors(strainCtr,:))
        gr = cat(2,gr{:});
        boundedline(distBins(2:end)-distBinwidth/2,nanmean(gr,2),...
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
        figurename = ['figures/speedvsneighbrdistance_' wormnum{1}];
        exportfig(speedFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
%         dircorrFig.Children.YLim = [-1 1];
        dircorrFig.Children.XLim = [0 2000];
        set(dircorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
        ylabel(dircorrFig.Children,'orientational correlation')
        xlabel(dircorrFig.Children,'distance between pair (μm)')
        legend(dircorrFig.Children,lineHandles,strains)
        figurename = ['figures/dircrosscorr_' wormnum{1}];
        exportfig(dircorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
%         velcorrFig.Children.YLim = [-1 1];
        velcorrFig.Children.XLim = [0 2000];
        set(velcorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
        ylabel(velcorrFig.Children,'velocity correlation')
        xlabel(velcorrFig.Children,'distance between pair (μm)')
        legend(velcorrFig.Children,lineHandles,strains)
        figurename = ['figures/velcrosscorr_' wormnum{1}];
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
        figurename = ['figures/radialdistributionfunction_' wormnum{1}];
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
