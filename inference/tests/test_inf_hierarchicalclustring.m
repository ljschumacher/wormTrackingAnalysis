% test inf_hierarchicalclustering
% test that hierarchical clustering calc behaves as expected for uniform random positions
% in both experiment and simulation modes

addpath('../component_functions')
L = 7.5;
numSamples = 2e4;
numWorms = 40;

%% test uniform random conditions
%% test for experimental conditions
% generate some uniform random coordinates
pix2mm = 0.0051;
x = L*rand(numSamples,1)/pix2mm;
y = L*rand(numSamples,1)/pix2mm;
frames = randi(numSamples/numWorms,numSamples,1);
branch_hist = inf_hierarchicalclustering({x,y,frames},'experiment',1);
plot(branch_hist)
hold on

%% test for simulation conditions
% generate some uniform random coordinates
x = reshape(L*rand(numSamples,1),numWorms,1,1,[]);
y = reshape(L*rand(numSamples,1),numWorms,1,1,[]);
branch_hist = inf_hierarchicalclustering(cat(3,x,y),'simulation-test',1);
plot(branch_hist)
hold on

%% test clustered conditions
L_clust = 1.25;
frac_in_cluster = 1/2;
% random location of cluster
x_cluster = (L - L_clust)*rand();
y_cluster = (L - L_clust)*rand();
%% test for experimental conditions
% fraction of the worms are in the cluster region
x_clustw = (x_cluster + L_clust*rand(numSamples*frac_in_cluster,1))/pix2mm;
y_clustw = (y_cluster + L_clust*rand(numSamples*frac_in_cluster,1))/pix2mm;
x = [x_clustw; L*rand(numSamples*(1-frac_in_cluster),1)/pix2mm];
y = [y_clustw; L*rand(numSamples*(1-frac_in_cluster),1)/pix2mm];
frames = randi(numSamples/numWorms,numSamples,1);
branch_hist = inf_hierarchicalclustering({x,y,frames},'experiment',1);
plot(branch_hist)
hold on

%% test for simulation conditions
% fraction of the worms are in the cluster region
x_clustw = x_cluster + L_clust*rand(numSamples*frac_in_cluster,1);
y_clustw = y_cluster + L_clust*rand(numSamples*frac_in_cluster,1);
x = reshape([x_clustw; L*rand(numSamples*(1-frac_in_cluster),1)],numWorms,1,1,[]);
y = reshape([y_clustw; L*rand(numSamples*(1-frac_in_cluster),1)],numWorms,1,1,[]);
branch_hist = inf_hierarchicalclustering(cat(3,x,y),'simulation-test',1);
plot(branch_hist)
hold on
