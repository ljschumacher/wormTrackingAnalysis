% script for loading HDF5 trajectory files from a list, counting the number of worms
% i.e. tracks per frame, and plot a histogram against frame number. This
% script filters data based on intensity thresholding to remove larvae, and writes the
% logical index (as int32) by intensity thresholding into the skeleton hdf5 file under
% "/intensity_threshold_logical_index".

% set intensity threshold. Set to 100 for single worm, 50 for 40 worm, and
% 40 for high density movies. 
IntensityThres = 50;

% retrieve the list of files to plot, one line at a time
fid = fopen('list2plot.txt');
sfpath = fgetl(fid);


% set counter for assigning line color
ii = 1;

while ischar(sfpath)
    disp(sfpath)
    
    % load current trajectory data
    trajData = h5read(sfpath,'/trajectories_data');
    
    % remove data by intensity threshold
    BlobFeats = h5read(sfpath,'/blob_features');
    ValidWormIndex = int32(BlobFeats.intensity_mean > IntensityThres);
    %ValidWormIndex = h5read(sfpath,'/ValidWormIndex_IntensityThreshold');
    Frames = trajData.frame_number;
    Frames = Frames .* ValidWormIndex;
    lowIntIndices = find(Frames == 0);
    Frames(lowIntIndices) = [];
    %% write intensity logical index into the hdf5
    %fidd = H5F.open(sfpath,'H5F_ACC_RDWR','H5P_DEFAULT');
    %if H5L.exists(fidd,'/ValidWormIndex_IntensityThreshold','H5P_DEFAULT')
    %   H5L.delete(fidd,'/ValidWormIndex_IntensityThreshold','H5P_DEFAULT');
    %end
    %H5F.close(fidd)
    %h5create(sfpath,'/ValidWormIndex_IntensityThreshold', size(ValidWormIndex),'Datatype','int32')
    %h5write(sfpath,'/ValidWormIndex_IntensityThreshold', ValidWormIndex)
  
    %%
    % plot histogram
    plotcolor = colorcube(15);
    histogram(Frames,'BinWidth',9,'DisplayStyle','stairs','EdgeColor',plotcolor(ii,:),'Normalization','countdensity')  
    xlabel('frame number','FontSize',20)
    ylabel('number of tracked objects','FontSize',20)
    set(gca,'FontSize',15)
    hold on
    
    % go to the next line/file
    sfpath = fgetl(fid);
    ii = ii+1;
end
legend('1','2','3','4','5','6','7','8','9','10','11')
fclose(fid);
fprintf 'Done plotting\n'
