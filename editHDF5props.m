

clear
close all

strains = {'npr1','N2','HA'};
wormnums = {'1W','40','HD'};
intensityThresholds = [100, 50, 40];

for strainCtr = 1:length(strains)
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_list.txt']);
        numFiles = length(filenames);
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
%             try
                frameRate = h5readatt(filename,'/plate_worms','expected_fps');
%             catch
% %                 frameRate = input(['Please enter /expected_fps for ' filename]);
%                 frameRate = h5readatt(filename,'/plate_worms','/expected_fps');
%                 h5writeatt(filename,'/plate_worms','expected_fps',frameRate)
%             end
            if ~ismember(frameRate,[3.0, 9.0])
                frameRate = input(['Please enter expected_fps for ' filename...
                    ', currently set to ' num2str(frameRate) ':']);
                h5writeatt(filename,'/plate_worms','expected_fps',frameRate)
            end
        end
    end
   
end