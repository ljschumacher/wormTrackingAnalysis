%% generate tif images for hand labeling
    % get path to the MaskedVideo file
    maskedVideoFileName = '/Volumes/behavgenom_archive$/Serena/bigpatch/MaskedVideos/bigpatch_1.1_comp1_DA609_N2/bigpatch_1.1_comp1_DA609_N2_Ch1_04062018_162821.hdf5';
    % read full image from the MaskedVideo
    fullData = h5read(maskedVideoFileName,'/full_data');
    firstFullImage = fullData(:,:,1);
    % save first image
    splitMaskedVideoFileName = strsplit(maskedVideoFileName,'/');
    imageFileName1 = splitMaskedVideoFileName{end-1};
    imageFileName2 = splitMaskedVideoFileName{end};
    imageFileName2 = strrep(imageFileName2,'.hdf5','.jpg');
    imageFileName = ['/Volumes/behavgenom_archive$/Serena/AggregationScreening/manualFoodContourImages/' imageFileName1 '__' imageFileName2];
    imwrite(firstFullImage,imageFileName);
%% hand label food contour using VGG annotator (http://www.robots.ox.ac.uk/~vgg/software/via/via.html) and save annotations