% calculate speed vs neighbour distance, directional correlation, and
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
% which are 1.57*iqr/sqrt(n)
distBinwidth = 5; % in units of micrometers
maxDist = 2000;
distBins = 0:distBinwidth:maxDist;
speedxticks = 0:250:2000;
dircorrxticks = 0:50:250;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'npr1','N2'};
wormnums = {'40','HD'};
intensityThresholds = [60, 40];
maxBlobSize = 1e4;
plotDiagnostics = false;
visitfreqFig = figure; hold on
for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        speedFig = figure; hold on
        dircorrFig = figure; hold on
        poscorrFig = figure; hold on
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
        numFiles = length(filenames);
        if wormnum == '40'
            visitfreq = cell(numFiles,1);
        end
        speeds = cell(numFiles,1);
        dxcorr = cell(numFiles,1); % for calculating directional cross-correlation
        pairdist = cell(numFiles,1);
        mindist= cell(numFiles,1);
        gr =cell(numFiles,1);
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate/10);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
            %% filter worms
            if plotDiagnostics
                plotIntensitySizeFilter(blobFeats,pixelsize,...
                    intensityThresholds(numCtr),maxBlobSize,...
                    [wormnum ' ' strains{strainCtr} ' ' strrep(filename(end-38:end-23),'/','')])
            end
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            %% calculate stats
            speeds{fileCtr} = cell(numFrames,1);
            dxcorr{fileCtr} = cell(numFrames,1); % for calculating directional cross-correlation
            pairdist{fileCtr} = cell(numFrames,1);
            mindist{fileCtr}= cell(numFrames,1);
            gr{fileCtr} = NaN(length(distBins) - 1,numFrames);
            for frameCtr = 1:numFrames
                frame = framesAnalyzed(frameCtr);
                [x, y, u, v] = calculateWormSpeeds(trajData, frame);
                if numel(x)>1 % need at least two worms in frame
                    speeds{fileCtr}{frameCtr} = sqrt(u.^2+v.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
                    dxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(u,v,true); % directional correlation
                    pairdist{fileCtr}{frameCtr} = pdist([x y]).*pixelsize; % distance between all pairs, in micrometer
                    gr{fileCtr}(:,frameCtr) = histcounts(pairdist{fileCtr}{frameCtr},distBins,'Normalization','probability'); % radial distribution function
                    gr{fileCtr}(:,frameCtr) = gr{fileCtr}(:,frameCtr)'.*maxDist^2./(2*distBins(2:end)*distBinwidth); % normalization
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
            pairdist{fileCtr} = horzcat(pairdist{fileCtr}{:});
            mindist{fileCtr} = horzcat(mindist{fileCtr}{:});
            % heat map of sites visited - this only makes sense for 40 worm
            % dataset where we don't move the camera
            if strcmp(wormnum,'40')&& plotDiagnostics
                siteVisitFig = figure;
                h=histogram2(trajData.coord_x*pixelsize/1000,trajData.coord_y*pixelsize/1000,...
                    'DisplayStyle','tile','EdgeColor','none','Normalization','probability');
                visitfreq{fileCtr} = h.Values(:);
                cb = colorbar; cb.Label.String = '# visited';
                axis equal
                xlabel('x (mm)'), ylabel('y (mm)')
                title([strains{strainCtr} ' ' strrep(filename(end-38:end-23),'/','')])
                set(siteVisitFig,'PaperUnits','centimeters')
                figurename = ['figures/individualRecordings/' strains{strainCtr} '_' strrep(strrep(filename(end-38:end-23),' ',''),'/','') '_sitesVisited'];
                exportfig(siteVisitFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
        [s_med,s_mad] = grpstats(vertcat(speeds{:}),quant(horzcat(mindist{:}),distBinwidth),...
            {@median,mad1});
        [c_med,c_mad] = grpstats(horzcat(dxcorr{:}),quant(horzcat(pairdist{:}),distBinwidth),...
            {@median,mad1});
        %% plot data
        bins = (0:numel(s_med)-1).*distBinwidth;
        boundedline(bins,smooth(s_med),[smooth(s_mad), smooth(s_mad)],...
            'alpha',speedFig.Children)
        bins = (0:numel(c_med)-1).*distBinwidth;
        boundedline(bins,smooth(c_med),[smooth(c_mad), smooth(c_mad)],...
            'alpha',dircorrFig.Children)
        gr = cat(2,gr{:});
        boundedline(distBins(2:end)-distBinwidth/2,mean(gr,2),[mad(gr,0,2) mad(gr,0,2)],...
            'alpha',poscorrFig.Children)
        if  strcmp(wormnum,'40')&& plotDiagnostics
            histogram(visitfreqFig.Children,vertcat(visitfreq{:}),'DisplayStyle','stairs','Normalization','probability')
        end
        %% format and export figures
        for figHandle = [speedFig, dircorrFig, poscorrFig] % common formating for both figures
            title(figHandle.Children,strains{strainCtr},'FontWeight','normal');
            set(figHandle,'PaperUnits','centimeters')
        end
        %
        speedFig.Children.YLim = [0 400];
        speedFig.Children.XLim = [0 2000];
        speedFig.Children.Box = 'on';
        speedFig.Children.XDir = 'reverse';
        ylabel(speedFig.Children,'speed (\mum/s)')
        xlabel(speedFig.Children,'distance to nearest neighbour (\mum)')
        figurename = ['figures/' strains{strainCtr} '_' wormnum '_speedvsneighbourdistance'];
        exportfig(speedFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
        dircorrFig.Children.YLim = [-1 1];
        dircorrFig.Children.XLim = [0 250];
        set(dircorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
        ylabel(dircorrFig.Children,'directional cross-correlation')
        xlabel(dircorrFig.Children,'distance between pair (\mum)')
        figurename = ['figures/' strains{strainCtr} '_' wormnum '_dircrosscorr'];
        exportfig(dircorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
        poscorrFig.Children.YLim = [0 10];
        poscorrFig.Children.XLim = [0 2000];
        poscorrFig.Children.Box = 'on';
        ylabel(poscorrFig.Children,'radial distribution function g(r)')
        xlabel(poscorrFig.Children,'distance r (\mum)')
        figurename = ['figures/' strains{strainCtr} '_' wormnum '_radialdistributionfunction'];
        exportfig(poscorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
    if  strcmp(wormnum,'40')&& plotDiagnostics
        visitfreqFig.Children.XScale = 'log';
        visitfreqFig.Children.YScale = 'log';
        visitfreqFig.Children.XLim = [4e-5 1e-1];
        xlabel(visitfreqFig.Children,'site visit frequency, f')
        ylabel(visitfreqFig.Children,'probability p(f)')
        legend(visitfreqFig.Children,strains)
        figurename = ['figures/' wormnum '_visitfreq'];
        exportfig(visitfreqFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end