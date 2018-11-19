% A. Script visualises sweeping using pixel data from selected frames, using up to 5 hours of long sweeping brightfield recordings (Leah's 1 patch dataset)
% step 1: read image frame
% step 2: generate binary image based on intensity threshold to pick out worm pixels from background
% step 3: dilate image to "connect" loose pharynxes and apply area thresholding to binary image to pick up clusters
% step 4: draw clusters over time (optional: plot centroid)
% step 5: plot food contour on top of the image (optional)
% B,C. Script calculates cluster centroid speed and plots them as timeseries and box plot.
% D. Script plots MSD of cluster centroids
% E. Script generates a down-sampled avi video, stringing together segments over multiple 1-hour recordings


close all
clear

bigpatch = true;
strains = {'npr1','N2'};
sampleEveryNSec = 120;  % in seconds
blobAreaThreshold = 3000; % single worm area ~ 500
plotVisualisation = true;
saveCentroidValues = false;
plotCentroidSpeeds = false;
plotMeanSquaredDisplacement = false;
makeDownSampledVideo = false;

if ~bigpatch
    maxSeg = 5;
else
    maxSeg = 20;
end

if plotVisualisation
    plotCentroid = false;
    plotFoodContour = true;
end

if plotCentroidSpeeds || plotMeanSquaredDisplacement
    smoothWindow = 20; % number of sampled frames for smoothing i.e. smoothWindow = 20x sampleEveryNSec = 30 means smoothing over 600 seconds
    if ~bigpatch
        plotFileList = [1,2,3,4,6]; % the replicates to plot for overall cluster speed - ignore rep 5 with two converging clusters
    end
end

if plotMeanSquaredDisplacement
    initialTimeStepToUse = 60; % use file matching this timeStep to start with
end

pixelsize = 10; % microns per pixel. can read by pixelsize = double(h5readatt(skelFilename,'/trajectories_data','microns_per_pixel'));
frameRate = 25; % frames per second.

exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

exportOptions2 = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',25,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1); % for individual rep timeseries only

addpath('auxiliary/')

for strainCtr = 1%:length(strains)
    if ~bigpatch
        [annotationNum,annotationFilenames,~] = xlsread('datalists/BFLongSweeping.xlsx',strainCtr,'A1:E30','basic');
    else
        [annotationNum,annotationFilenames,~] = xlsread('datalists/BFBigPatchLongSweeping.xlsx',strainCtr,'A1:E200','basic');
    end
    % xy coordinates and radius of food contour obtained by hand annotation using VGG
    if ~bigpatch
        if strcmp(strains{strainCtr},'npr1')
            foodCtnCoords_xyr = [1203,914,379;1057,867,348;997,810,328;1007,790,334;988,711,335;1006,683,327];
        elseif strcmp(strains{strainCtr},'N2')
            foodCtnCoords_xyr = [1055,800,380;1194,754,356;714,783,322;1328,881,338;1022,678,328;905,786,330];
        end
    else
        if strcmp(strains{strainCtr},'npr1')
            foodCtnCoords_x{1} = [427,417,433,482,526,620,706,809,863,950,1024,1107,1226,1281,1393,1502,1554,1698,1733,1746,1753,1762,1766,1746,1717,1653,1534,1428,1313,1207,1124,1043,947,841,770,703,607,559,510,472,437,427];
            foodCtnCoords_x{2} = [300,312,330,375,414,468,543,681,750,894,960,1043,1136,1208,1301,1445,1493,1595,1634,1637,1640,1613,1562,1445,1382,1292,1169,1040,936,837,741,651,567,507,363,330,306,300];
            foodCtnCoords_x{3} = [222,264,330,393,498,591,690,819,918,1040,1145,1241,1307,1409,1523,1544,1565,1580,1577,1475,1403,1298,1166,1031,915,786,657,528,435,303,252,228,219,216,222];
            foodCtnCoords_x{4} = [361,363,396,445,493,557,620,697,776,876,953,1018,1112,1199,1264,1336,1411,1476,1647,1687,1710,1712,1707,1652,1585,1498,1401,1294,1180,1045,963,824,714,642,562,483,445,408,383,366,361];
            foodCtnCoords_x{5} = [276,301,351,413,470,570,669,766,911,995,1105,1234,1341,1438,1510,1548,1565,1573,1573,1555,1493,1421,1304,1194,1045,943,814,702,605,490,413,346,304,276,259,256,276];
            foodCtnCoords_x{6} = [162,179,226,279,326,396,468,555,642,724,844,970,1063,1197,1292,1366,1428,1471,1520,1528,1528,1510,1458,1391,1316,1234,1155,1030,928,814,677,565,435,353,284,236,192,174,162];
            foodCtnCoords_x{7} = [421,458,498,570,664,729,816,923,1058,1185,1366,1436,1583,1615,1652,1652,1632,1553,1406,1202,1127,1045,913,791,659,577,518,473,428,403,388,378,378,421];
            foodCtnCoords_y{1} = [995,1120,1255,1384,1454,1566,1611,1650,1676,1695,1701,1705,1692,1688,1650,1579,1522,1287,1204,1130,1046,960,883,796,713,613,510,446,392,385,376,382,398,437,472,517,610,674,751,835,924,995];
            foodCtnCoords_y{2} = [1020,1139,1274,1358,1427,1487,1547,1640,1670,1703,1709,1703,1679,1667,1634,1514,1463,1289,1184,1049,945,822,729,576,525,471,417,387,378,396,411,447,510,576,774,855,936,1020];
            foodCtnCoords_y{3} = [1088,1292,1391,1481,1574,1628,1670,1697,1706,1691,1658,1610,1559,1463,1271,1193,1109,996,852,621,537,456,378,351,339,363,387,447,525,696,828,918,987,1067,1088];
            foodCtnCoords_y{4} = [983,1082,1209,1309,1396,1486,1525,1580,1613,1647,1662,1670,1667,1652,1640,1610,1560,1508,1274,1180,1080,951,834,704,590,500,418,368,328,321,328,358,403,453,525,612,692,771,854,928,983];
            foodCtnCoords_y{5} = [1063,1172,1269,1351,1421,1493,1553,1583,1605,1603,1593,1543,1491,1379,1296,1187,1080,970,871,747,577,470,383,323,279,264,276,311,346,426,508,617,704,806,891,995,1063];
            foodCtnCoords_y{6} = [1005,1152,1294,1371,1443,1506,1563,1610,1635,1655,1680,1680,1667,1627,1555,1476,1406,1309,1197,1075,975,863,744,647,575,515,458,418,401,393,401,433,503,570,654,724,836,928,1005];
            foodCtnCoords_y{7} = [1292,1394,1461,1565,1625,1660,1702,1729,1734,1725,1630,1575,1371,1247,1102,946,834,667,508,388,363,341,351,393,453,533,600,674,769,844,958,1085,1185,1292];
        elseif strcmp(strains{strainCtr},'N2')
            foodCtnCoords_x{1} = [315,312,315,327,354,387,429,486,585,651,756,828,921,1043,1142,1229,1367,1499,1589,1652,1643,1580,1478,1367,1274,1142,918,786,669,597,513,462,387,333,315];
            foodCtnCoords_x{2} = [441,465,492,531,603,666,762,858,936,993,1073,1181,1232,1328,1388,1502,1565,1712,1775,1826,1853,1862,1862,1856,1799,1655,1517,1445,1361,1229,1088,939,789,690,624,561,471,444,438,441];
            foodCtnCoords_x{3} = [381,369,375,399,441,483,624,684,774,861,999,1085,1193,1310,1409,1481,1562,1637,1664,1676,1694,1688,1655,1601,1538,1466,1385,1298,1202,1073,966,849,729,618,522,468,417,381];
            foodCtnCoords_x{4} = [440,445,463,500,540,607,699,769,901,1015,1130,1264,1354,1426,1496,1550,1640,1685,1739,1779,1819,1829,1824,1802,1762,1712,1615,1540,1428,1294,1162,1013,873,776,682,600,530,485,465,445,440];
            foodCtnCoords_x{5} = [500,500,508,547,602,657,732,834,928,1063,1142,1217,1316,1401,1493,1578,1645,1692,1727,1737,1737,1707,1647,1555,1456,1356,1217,1078,973,898,809,719,647,585,528,505,500];
            foodCtnCoords_x{6} = [244,241,271,314,371,411,485,572,642,722,841,968,1063,1167,1259,1326,1406,1468,1535,1583,1613,1620,1610,1590,1540,1473,1336,1219,1140,1043,953,886,766,657,565,488,431,378,338,304,269,254,244];
            foodCtnCoords_x{7} = [431,428,460,565,625,714,839,916,1048,1150,1257,1349,1466,1540,1615,1670,1729,1754,1754,1747,1722,1680,1608,1533,1379,1279,1170,1018,868,754,647,577,515,480,450,431];
            foodCtnCoords_y{1} = [1023,1172,1271,1364,1445,1526,1577,1643,1727,1763,1802,1814,1820,1802,1787,1754,1667,1511,1319,1097,1014,753,609,513,462,420,402,438,483,534,621,681,792,933,1023];
            foodCtnCoords_y{2} = [1076,1184,1277,1349,1433,1499,1574,1619,1643,1655,1661,1661,1658,1625,1607,1544,1493,1346,1259,1130,1040,957,855,795,618,471,390,366,345,339,330,357,420,459,513,633,816,912,1008,1076];
            foodCtnCoords_y{3} = [996,1109,1259,1352,1451,1520,1652,1691,1730,1763,1796,1793,1787,1736,1670,1628,1544,1418,1346,1247,1127,1034,885,768,699,624,573,531,507,480,486,504,558,618,723,816,894,996];
            foodCtnCoords_y{4} = [1172,1254,1346,1448,1508,1588,1670,1715,1772,1807,1814,1799,1789,1759,1722,1687,1617,1553,1471,1391,1274,1165,1080,980,878,784,679,622,557,518,505,508,545,585,647,739,821,928,1028,1112,1172];
            foodCtnCoords_y{5} = [1048,1170,1282,1379,1478,1553,1625,1697,1737,1767,1764,1754,1727,1690,1642,1563,1468,1379,1262,1137,1025,876,756,642,557,505,465,458,478,498,535,607,674,744,863,980,1048];
            foodCtnCoords_y{6} = [1090,1187,1341,1448,1523,1578,1655,1712,1747,1764,1784,1787,1782,1752,1707,1675,1608,1538,1446,1324,1212,1078,961,868,774,659,545,483,455,438,433,435,448,483,550,607,667,749,816,893,993,1048,1090];
            foodCtnCoords_y{7} = [1115,1254,1408,1605,1662,1737,1804,1822,1834,1836,1819,1797,1739,1670,1595,1515,1389,1237,1130,1023,918,834,739,652,570,542,518,515,557,612,689,761,871,941,1023,1115];
        end
    end
    
    % go through each recording replicate (6 in total)
    for fileCtr = 7%1:6
        if plotVisualisation
            clusterVisFig = figure; hold on
        end
        if makeDownSampledVideo
            if ~bigpatch
                video = VideoWriter([strains{strainCtr} '_rep' num2str(fileCtr) '_timeStep' num2str(sampleEveryNSec) '.avi']); % create the video object
            else
                video = VideoWriter([strains{strainCtr} '_bigpatch_rep' num2str(fileCtr) '_timeStep' num2str(sampleEveryNSec) '.avi']); % create the video object
            end
            video.FrameRate = 15;% set the frame rate
            open(video); % open the file for writing
        end
        totalFrames = 0;
        totalSegs = 0;
        for segCtr = 1:maxSeg % go through each hour of the recording replicate (5 hours maximum)
            fileIdx = find(annotationNum(:,1) == fileCtr & annotationNum(:,2) == segCtr);
            firstFrame{segCtr} = annotationNum(fileIdx,4)+1; % +1 to adjust for python 0 indexing
            lastFrame{segCtr} = annotationNum(fileIdx,5)+1;
            filename{segCtr} = annotationFilenames{fileIdx};
            if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
                totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
                totalSegs = totalSegs+1;
            end
        end
        totalSampleFrames = ceil(totalFrames/sampleEveryNSec/25);
        % initialise
        plotColors = parula(totalSampleFrames);
        clusterCentroidCoords{fileCtr} = NaN(2,totalSampleFrames);
        clusterSolidity{fileCtr} = NaN(1,totalSampleFrames);
        clusterArea{fileCtr} = NaN(1,totalSampleFrames);
        cumFrame = 0; % keep track of cumulative frames across replicate segments
        leftoverFrames = 0; % keep track of leftover frames at the end of one segment that combines with the start of the next segment
        
        for segCtr = 1:totalSegs
            % load data
            %skelFilename = strrep(strrep(filename{segCtr},'MaskedVideos','Results'),'.hdf5','_skeletons.hdf5');
            %frameRate = double(h5readatt(skelFilename,'/plate_worms','expected_fps'));
            %trajData = h5read(skelFilename,'/trajectories_data');
            fileInfo = h5info(filename{segCtr});
            dims = fileInfo.Datasets(2).Dataspace.Size;
            if leftoverFrames>0
                assert(firstFrame{segCtr} ==1); % if there are leftover frames from the previous segment, then this segment must start from the very first frame
                firstFrame{segCtr} = sampleEveryNSec*frameRate-leftoverFrames+firstFrame{segCtr};
            end
            movieFrames = firstFrame{segCtr}:sampleEveryNSec*frameRate:...
                floor((lastFrame{segCtr}-firstFrame{segCtr}+1)/sampleEveryNSec/frameRate)*sampleEveryNSec*frameRate+firstFrame{segCtr};
            for frameCtr = 1:numel(movieFrames)
                imageFrame = h5read(filename{segCtr},'/mask',[1,1,movieFrames(frameCtr)],[dims(1),dims(2),1]);
                if makeDownSampledVideo
                    timeStamp = round((cumFrame+frameCtr-1)*sampleEveryNSec/60); % timestamp in minutes
                    imageFrame = AddTextToImage(imageFrame,['t = ' num2str(timeStamp) ' min'],[100 100],[255,255,255],'Arial',50);
                    writeVideo(video,imageFrame); %write the image to file
                end
                % generate binary segmentation based on black/white contrast
                binaryImage = imageFrame>0 & imageFrame<70;
                binaryImage = imfill(binaryImage, 'holes');
                % filter by blob size
                blobMeasurements = regionprops(binaryImage, 'Area','Centroid','Solidity');
                blobCentroidsCoords = reshape([blobMeasurements.Centroid],[2, numel([blobMeasurements.Centroid])/2]);
                blobLogInd = [blobMeasurements.Area] > blobAreaThreshold; % apply blob area threshold values
                blobLogInd = blobLogInd & blobCentroidsCoords(2,:) > 250; % get rid of the annoying box at the edge
                % restrict to blobs near the food patch centre (within 500 pixels or 5 mm)
                if ~bigpatch
                    blobLogInd = blobLogInd & blobCentroidsCoords(1,:)<foodCtnCoords_xyr(fileCtr,1)+500 & blobCentroidsCoords(1,:)>foodCtnCoords_xyr(fileCtr,1)-500;
                    blobLogInd = blobLogInd & blobCentroidsCoords(2,:)<foodCtnCoords_xyr(fileCtr,2)+500 & blobCentroidsCoords(2,:)>foodCtnCoords_xyr(fileCtr,2)-500;
                end
                blobBoundaries = bwboundaries(binaryImage,8,'noholes');
                if plotVisualisation
                    % plot individual blob boundaries that meet area threshold requirements
                    set(0,'CurrentFigure',clusterVisFig)
                    for blobCtr = 1:numel(blobLogInd)
                        if blobLogInd(blobCtr)
                            fill(blobBoundaries{blobCtr}(:,1)*pixelsize/1000,blobBoundaries{blobCtr}(:,2)*pixelsize/1000,plotColors(cumFrame+frameCtr,:),'edgecolor','none')
                            alpha 0.5
                        end
                    end
                end
                % get centroids and solidity
                if nnz(blobLogInd)>0
                    [~,maxAreaIdx] = max([blobMeasurements.Area].*blobLogInd);
                    clusterCentroidCoords{fileCtr}(1,cumFrame+frameCtr) = blobCentroidsCoords(1,maxAreaIdx);
                    clusterCentroidCoords{fileCtr}(2,cumFrame+frameCtr) = blobCentroidsCoords(2,maxAreaIdx);
                    area = [blobMeasurements.Area];
                    solidity = [blobMeasurements.Solidity];
                    clusterSolidity{fileCtr}(cumFrame+frameCtr) = solidity(maxAreaIdx);
                    clusterArea{fileCtr}(cumFrame+frameCtr) = area(maxAreaIdx);
                    if plotVisualisation & plotCentroid
                        set(0,'CurrentFigure',clusterVisFig)
                        plot(clusterCentroidCoords{fileCtr}(2,cumFrame+frameCtr)*pixelsize/1000,clusterCentroidCoords{fileCtr}(1,cumFrame+frameCtr)*pixelsize/1000,'k--x')
                    end
                end
            end
            cumFrame = cumFrame+numel(movieFrames);
            leftoverFrames = lastFrame{segCtr} - max(movieFrames);
        end
        
        if plotVisualisation
            set(0,'CurrentFigure',clusterVisFig)
            axis equal
            colorbar
            caxis([0 ceil(totalFrames/25/60)])
            cb = colorbar; cb.Label.String = 'minutes';
            if ~bigpatch
                xmax = round(foodCtnCoords_xyr(fileCtr,2)*pixelsize/1000+5);
                xmin = round(foodCtnCoords_xyr(fileCtr,2)*pixelsize/1000-5);
                ymax = round(foodCtnCoords_xyr(fileCtr,1)*pixelsize/1000+5);
                ymin = round(foodCtnCoords_xyr(fileCtr,1)*pixelsize/1000-5);
                xlim([xmin xmax])
                ylim([ymin ymax])
                xticks(xmin:2:xmax)
                yticks(ymin:2:ymax)
            else
                xlim([0 20]);
                ylim([0 20]);
                xticks([0:4:20])
                yticks([0:4:20])
            end
            xlabel('x (mm)')
            ylabel('y (mm)')
            if plotFoodContour
                if ~bigpatch
                    viscircles([foodCtnCoords_xyr(fileCtr,2),foodCtnCoords_xyr(fileCtr,1)]*pixelsize/1000,foodCtnCoords_xyr(fileCtr,3)*pixelsize/1000,'Color','k','LineStyle','--','LineWidth',1);
                else
                    plot(foodCtnCoords_x{fileCtr}*pixelsize/1000,foodCtnCoords_y{fileCtr}*pixelsize/1000,'Color','k','LineStyle','--','LineWidth',1);
                end
            end
            % export figure
            if ~bigpatch
                if plotCentroid
                    if plotFoodContour
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroidFood_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
                    else
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroid_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
                    end
                elseif plotFoodContour
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelFood_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
                else
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixel_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
                end
            else
                if plotCentroid
                    if plotFoodContour
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroidFood_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFBP'];
                    else
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroid_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFBP'];
                    end
                elseif plotFoodContour
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelFood_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFBP'];
                else
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixel_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFBP'];
                end
            end
            exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
        end
        
        % calculate centroid speed (in microns per minute)
        clusterCentroidSpeed{fileCtr} = NaN(1,totalSampleFrames);
        for frameCtr = 1:totalSampleFrames-10
            if ~isnan(clusterCentroidCoords{fileCtr}(1,frameCtr))
                for stepCtr = 1:10 % in case the next sample frame has no cluster, go up to 10 time steps away
                    if ~isnan(clusterCentroidCoords{fileCtr}(1,frameCtr+stepCtr))
                        break
                    end
                end
                clusterCentroidSpeed{fileCtr}(frameCtr) = sqrt((clusterCentroidCoords{fileCtr}(1,frameCtr+stepCtr)-clusterCentroidCoords{fileCtr}(1,frameCtr))^2 +...
                    (clusterCentroidCoords{fileCtr}(2,frameCtr+stepCtr)-clusterCentroidCoords{fileCtr}(2,frameCtr))^2)...
                    /stepCtr*pixelsize*60/sampleEveryNSec;
                if clusterCentroidSpeed{fileCtr}(frameCtr)>500
                    clusterCentroidSpeed{fileCtr}(frameCtr) = NaN;
                end
            end
        end
        if makeDownSampledVideo
            close(video); %close the file
        end
    end
    
    % save cluster centroid coordinates and speeds
    
    if saveCentroidValues
        if ~bigpatch
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidSpeed');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidCoords_dataBF_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidCoords');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSolidity_dataBF_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterSolidity');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidArea_dataBF_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterArea');
        else
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_dataBFBP_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidSpeed');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidCoords_dataBFBP_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidCoords');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSolidity_dataBFBP_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterSolidity');
            save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidArea_dataBFBP_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterArea');
        end
    end
end

%% plot median speeds for different replicates (npr-1 only)
if plotCentroidSpeeds
    load(['/Users/sding/Documents/trackingAnalysis/figures/sweeping/npr1_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '.mat'])
    smoothSpeeds = NaN(numel(plotFileList),maxSeg*3600/sampleEveryNSec); % initialise
    recordingColors = distinguishable_colors(numel(plotFileList));
    legends = cell(1,numel(plotFileList));
    
    clusterCentroidSpeedFig = figure; hold on % show each replicate individually
    poolRepFig = figure; % shaded error bars showing average across each specified replicate
    smoothedBoxPlotFig = figure; % show each replicate as a box plot
    unsmoothedBlotPlotFig = figure; % show each replicate as a box plot
    speedVSolidityFig = figure; % plots speed against solidity
    
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        totalFrames = 0;
        totalSegs = 0;
        for segCtr = 1:5 % go through each hour of the recording replicate (5 hours maximum)
            recIdx = find(annotationNum(:,1) == fileCtr & annotationNum(:,2) == segCtr);
            firstFrame{segCtr} = annotationNum(recIdx,4)+1; % +1 to adjust for python 0 indexing
            lastFrame{segCtr} = annotationNum(recIdx,5)+1;
            filename{segCtr} = annotationFilenames{recIdx};
            if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
                totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
                totalSegs = totalSegs+1;
            end
        end
        recordingsPlotX = 1:sampleEveryNSec/60:ceil(totalFrames/25/60);
        set(0,'CurrentFigure',clusterCentroidSpeedFig)
        if numel(recordingsPlotX) == numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) < numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileIdx}(1:numel(recordingsPlotX)),'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) > numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX(1:numel(clusterCentroidSpeed{fileIdx})),smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        end
        legends{fileCtr} = num2str(fileIdx);
        
        if fileIdx ==2
            speedY = smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan');
            smoothSpeeds(fileCtr,1:numel(speedY)-20)=speedY(21:end);
        else
            smoothSpeeds(fileCtr,1:numel(clusterCentroidSpeed{fileIdx}))=...
                smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan');
        end
        %figure; plot(smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),clusterSolidity{fileIdx},'.')
    end
    
    set(0,'CurrentFigure',clusterCentroidSpeedFig)
    xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
    legend(legends,'Location','eastoutside')
    xlim([0 300])
    set(gca,'Xtick',[0:50:300])
    if ~bigpatch
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_individualReps'];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFBP_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_individualReps'];
    end
    exportfig(clusterCentroidSpeedFig,[figurename '.eps'],exportOptions2)
    
    set(0,'CurrentFigure',poolRepFig)
    recordingsPlotX = 1:sampleEveryNSec/60:300;
    shadedErrorBar(recordingsPlotX,nanmedian(smoothSpeeds,1),nanstd(smoothSpeeds),'k');
    xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
    xlim([0 250])
    set(gca,'Xtick',[0:50:250])
    if ~bigpatch
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_pooledReps'];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFBP_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_pooledReps'];
    end
    exportfig(poolRepFig,[figurename '.eps'],exportOptions)
    
    set(0,'CurrentFigure',smoothedBoxPlotFig)
    boxGroups = [];
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        if sampleEveryNSec == 30
            boxGroups = [boxGroups fileIdx*ones(1,599)];
        elseif sampleEveryNSec ==120
            boxGroups = [boxGroups fileIdx*ones(1,150)];
        end
    end
    boxplot(smoothSpeeds(:),boxGroups(:))
    set(0,'CurrentFigure',smoothedBoxPlotFig)
    xlabel('Replicate'), ylabel('Cluster Speed (microns/min)')
    set(gca,'XTickLabel',legends)
    if ~bigpatch
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_smoothedBoxPlot'];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFBP_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_smoothedBoxPlot'];
    end
    exportfig(smoothedBoxPlotFig,[figurename '.eps'],exportOptions)
    
    set(0,'CurrentFigure',unsmoothedBlotPlotFig)
    unsmoothedSpeeds = [];
    unsmoothedGroups = [];
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        unsmoothedSpeeds = [unsmoothedSpeeds clusterCentroidSpeed{fileIdx}];
        unsmoothedGroups = [unsmoothedGroups fileIdx*ones(1,numel(clusterCentroidSpeed{fileIdx}))];
    end
    boxplot(unsmoothedSpeeds(:),unsmoothedGroups(:))
    xlabel('Replicate'), ylabel('Cluster Speed (microns/min)')
    if ~bigpatch
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBF_timeStep' num2str(sampleEveryNSec) '_unsmoothedBlotPlot'];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFBP_timeStep' num2str(sampleEveryNSec) '_unsmoothedBlotPlot'];
    end
    exportfig(unsmoothedBlotPlotFig,[figurename '.eps'],exportOptions)
    
end

if plotMeanSquaredDisplacement
    % read xy coordinates from specified initial timeStep file
    load(['/Users/sding/Documents/trackingAnalysis/figures/sweeping/npr1_clusterCentroidCoords_dataBF_timeStep' num2str(initialTimeStepToUse) '.mat'])
    % initialise
    msdFig = figure; hold on
    legends = cell(1,numel(plotFileList));
    
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        % smooth xy coordinates
        smoothClusterCentroidCoords = smoothdata(clusterCentroidCoords{fileIdx},1,'movmedian',smoothWindow,'omitnan');
        % interpolate over NaN xy coordinates
        xcoords = naninterp(smoothClusterCentroidCoords(1,:))';
        ycoords = naninterp(smoothClusterCentroidCoords(2,:))';
        data=[xcoords ycoords];
        % initialise
        nData = size(data,1); %number of data points
        numberOfDeltaT = floor(nData/4); %for MSD, dt should be up to 1/4 of number of data points (Saxton, M. J. (1997). Single-particle tracking: The distribution of diffusion coeffcients. Biophys. J. 72, 1744?1753.);
        msd = zeros(numberOfDeltaT,3); %We'll store [mean, std, n]
        % calculate msd for all deltaT's
        for dt = 1:numberOfDeltaT
            deltaCoords = data(1+dt:end,1:2) - data(1:end-dt,1:2);
            squaredDisplacement = sum(deltaCoords.^2,2); % dx^2+dy^2
            msd(dt,1) = mean(squaredDisplacement)*pixelsize/1e6; % average in cm
            msd(dt,2) = std(squaredDisplacement)*pixelsize/1e6; % std in cm
            msd(dt,3) = length(squaredDisplacement); % n
        end
        % plot MSD
        plot((1:size(msd,1))*initialTimeStepToUse/60,msd(:,1))
        % generate legend
        legends{fileCtr} = num2str(fileIdx);
    end
    % format MSD plot
    xlabel('tau (min)'); ylabel('MSD (cm)')
    legend(legends,'Location','eastoutside')
    if ~bigpatch
        figurename = ['figures/sweeping/npr1_clusterMSD_dataBF_initialTimeStep' num2str(initialTimeStepToUse)];
    else
        figurename = ['figures/sweeping/npr1_clusterMSD_dataBFBP_initialTimeStep' num2str(initialTimeStepToUse)];
    end
    exportfig(msdFig,[figurename '.eps'],exportOptions)
end