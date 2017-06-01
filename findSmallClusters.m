% Filter for small clusters (2-4 neighbrs) using distance to 10 nearest neighbrs written to red hdf5 files. 
% Returns matrices D with logical indices for small clusters, and E for nnz values of D. 
% Generates bar graphs to compare small cluster frequencies across strains and densities.

close all
clear

%% set parameters
strains = {'npr1','N2'};
wormnums = {'40','HD'};
minIntensities_r = [40, 80]; % script takes 40 for all movies but recordings 54 and 55, which takes 80 because the dynamic ranges were different for those movies
maxBlobSize_r = 250000;
minSkelLength_r = 850;
pixelsize = 100/19.5; % 100 microns is 19.5 pixels
loneClusterRadius = 2000;
inClusterRadius = 500;
minNumNeighbrs = [2,3,4];

%% go through different strains, densities, and movies
for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        % load files
        filenames_r = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames_r);
        for fileCtr = 1:numFiles
            filename_r = filenames_r{fileCtr};
            trajData_r = h5read(filename_r,'/trajectories_data');
            blobFeats_r = h5read(filename_r,'/blob_features');
            skelData_r = h5read(filename_r,'/skeleton');
            numCloseNeighbr_r = h5read(filename_r,'/num_close_neighbrs');
            neighbrDist_r = h5read(filename_r,'/neighbr_distances');
            % filter worms
            if isempty(find(filename_r == 54)) || isempty(find(filename_r == 55)) == (3>2)
                minIntensity = minIntensities_r(1);
            else
                minIntensity = minIntensities_r(2);
            end
            skelLengths = sum(sqrt(sum((diff(skelData_r,1,2)*pixelsize).^2)));
            trajData_r.filtered = (blobFeats_r.intensity_mean >= minIntensity)&...
                (blobFeats_r.area*pixelsize^2 <= maxBlobSize_r)&...
                logical(trajData_r.is_good_skel)&...
                logical(skelLengths(:)>minSkelLength_r);
            % filter for small clusters and write logical indices into D
            numNeighbrs = length(minNumNeighbrs);
            D = zeros(length(trajData_r.filtered),numNeighbrs);
            for neighbrCtr = 1:length(minNumNeighbrs);
                neighbrNum = minNumNeighbrs(neighbrCtr);
                D(:,neighbrCtr) = trajData_r.filtered&...
                    numCloseNeighbr_r== neighbrNum&...
                    neighbrDist_r(:,neighbrNum+1)>=loneClusterRadius;
            end
            % write number of clusters with 2-4 neighbors into E
            E(strainCtr,numCtr,fileCtr,1)=nnz(D(:,1));
            E(strainCtr,numCtr,fileCtr,2)=nnz(D(:,2));
            E(strainCtr,numCtr,fileCtr,3)=nnz(D(:,3));
        end
    end
end
%% plot graph
for neighbrCtr = 1:length(minNumNeighbrs);
    neighbrNum = minNumNeighbrs(neighbrCtr);
    figure;
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        for strainCtr = 1:length(strains)
            strain = strains{strainCtr};
            subplot(2,2,(numCtr-1)*length(strains)+strainCtr)
            bar(squeeze(E(numCtr,strainCtr,:,neighbrCtr)))
            title([strain ' ' wormnum],'FontWeight','normal')
        end
    end
    figName = strcat('smallCluster_',num2str(minNumNeighbrs(neighbrCtr)),'neighbrs_loneRadius',num2str(loneClusterRadius),'.fig');
    savefig(figName)
    %close all
end