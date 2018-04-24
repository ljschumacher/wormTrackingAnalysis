% quantify how much data has been tracked, and how much filtered out at
% varous stages
close all
clear

strains = {'npr1','N2'};
wormnums = {'1W','40','HD'};

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        figure, hold on
        % load red channel file list
        filenames = importdata(['../datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            % check how many red worms have been appropriately skeletonized
            numTracks = numel(unique(trajData.worm_index_joined));
            numFrames = numel(unique(trajData.frame_number));
            numBlobs = numel(trajData.has_skeleton);
            numSkels = nnz(trajData.has_skeleton);
            numSkelsGood = nnz(trajData.is_good_skel);
            plot([numBlobs, numSkels, numSkelsGood]./numFrames)
        end
        title([strain ' ' wormnum],'FontWeight','normal')
        xticks([1 2 3])
        xticklabels({'trajectories','skeletons','good skeletons'})
    end
end
tilefigs