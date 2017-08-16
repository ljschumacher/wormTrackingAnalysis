samplePathFig = figure; hold on
sampleAngSpeedFig = figure; hold on
samplePaths = importdata('headAngSpeedSample_leaveCluster2.mat');
%samplePaths = headAngSpeedSample_loneWorm2;
sampleAngSpeed = importdata('sHeadAngleSpeedSample_leaveCluster2.mat');
%sampleAngSpeed = sHeadAngleSpeedSample_loneWorm2;
samplePaths = samplePaths(~cellfun('isempty',samplePaths)); 
samplePaths = reshape(samplePaths,length(samplePaths)/2,2);
trajSamples = randi(size(samplePaths,1),[10,1]);
for sampleCtr = 1:length(trajSamples)
    trajSample = trajSamples(sampleCtr);
    sample_xcoords = samplePaths{trajSample,1};
    sample_ycoords = samplePaths{trajSample,2};
    samplem_xcoords = mean(samplePaths{trajSample,1},2);
    samplem_ycoords = mean(samplePaths{trajSample,2},2);
    samplem_xcoords = samplem_xcoords - samplem_xcoords(1);
    samplem_ycoords = samplem_ycoords - samplem_ycoords(1);
    sampleangspeed = sampleAngSpeed{trajSample};
    set(0,'CurrentFigure',samplePathFig)
    plot(samplem_xcoords,samplem_ycoords)
    set(0,'CurrentFigure',sampleAngSpeedFig)
    plot(sampleangspeed)
    waitforbuttonpress
end