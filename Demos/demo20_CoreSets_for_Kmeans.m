%{
Following review paper of:
"?Practical Coreset Constructions for Machine Learning"
by Bachem, Lucic, Krause 2017; ?http://arxiv.org/abs/1703.06476
%}

% Needs pdist2_faster, kmeansPlusPlus, hungarian, bestMap, mnist_data_all.mat
addpath ~/Repos/randomized-algorithm-class/Code/
addpath('~/Google Drive/TeachingDocs/APPM4720_5720_Spring19_Randomized/Code');
load mnist_data_all
percentCorrect = @(labels1,labels2) length(find(labels1==labels2))/length(labels1);
%%
p           = size(Train,2);
K           = 10; % ask for 10 labels

ALGO_NAMES  = {'Kmeans','Kmeans++','Kmeans-Coresets-uniform','Kmeans-Coresets'};
[TrainError,TestError,Timing]  = deal(zeros(length(ALGO_NAMES),1));

ALGO        = 1;
tic
[IDX_Train, ClusterCenters]    = kmeans( Train, K );
Timing(ALGO) = toc;
[~,IDX_Test   ]     = pdist2_faster(ClusterCenters,Test,'squaredeuclidean','smallest',1);
IDX_Train_permuted  = bestMap( Train_labels, IDX_Train );
IDX_Test_permuted   = bestMap( Test_labels, IDX_Test );
TrainError(ALGO) = percentCorrect(IDX_Train_permuted,Train_labels);
TestError(ALGO)  = percentCorrect(IDX_Test_permuted,Test_labels);
%% Use K-means++
ALGO        = 2;
tic
ClusterCenters = kmeansPlusPlus( Train, K );
Timing(ALGO) = toc;
[Dist_Kpp,IDX_Train_Kpp   ] = pdist2_faster(ClusterCenters,Train,'squaredeuclidean','smallest',1);
[~,IDX_Test_Kpp   ]  = pdist2_faster(ClusterCenters,Test,'squaredeuclidean','smallest',1);

IDX_Train_permuted  = bestMap( Train_labels, IDX_Train_Kpp );
IDX_Test_permuted   = bestMap( Test_labels, IDX_Test_Kpp );

TrainError(ALGO) = percentCorrect(IDX_Train_permuted,Train_labels);
TestError(ALGO)  = percentCorrect(IDX_Test_permuted,Test_labels);

%% Use K-means to get a core-set
N   = size( Train, 1 );
M   = round( N/100 ); % size of the core-set
% The naive/uniform coreset works fine if M is about N/10
%   but if we start sub-sampling further, e.g., N/100,
%   then the fancier weighted coreset starts to show improvement
%   (a worst-case improvement; on some random samples, it works fine)
alpha   = 16*(log(K)+2);
clusterAvgDist  = zeros(K,1);
c               = mean( Dist_Kpp );
weights         = zeros(N,1);
weights         = weights + alpha*Dist_Kpp'/c;
for k = 1:K
    ind         = find( IDX_Train_Kpp == k );
    clusterSize     = length( ind );
    ci          = mean( Dist_Kpp( ind ) );
    weights( ind ) = weights( ind ) + 2*alpha*ci/(c*clusterSize) + 4*N/clusterSize;
end
weights     = weights/sum(weights);
histogram( weights );

naive_coreset   = randsample( N, M ); % uniform weights
coreset         = randsample( N, M, true, weights );
% Now, to really do core-sets, we also need to update the weights for each
% entry that is sampled, e.g., before, it was implicitly 1/N
%   Now, it's now 1/M, but rather 1/(M*N*weights)
%   Not sure how to do that with Lloyd's algorithm much less
%   having to write our own kmeans script, so just ignore...
%% Now, re-run Kmeans on these sampled data

ALGO        = 3;
tic
[~, ClusterCenters]    = kmeans( Train(naive_coreset,:), K );
Timing(ALGO) = toc;
[~,IDX_Train  ]     = pdist2_faster(ClusterCenters,Train,'squaredeuclidean','smallest',1);
[~,IDX_Test   ]     = pdist2_faster(ClusterCenters,Test,'squaredeuclidean','smallest',1);
IDX_Train_permuted  = bestMap( Train_labels, IDX_Train );
IDX_Test_permuted   = bestMap( Test_labels, IDX_Test );
TrainError(ALGO)    = percentCorrect(IDX_Train_permuted,Train_labels);
TestError(ALGO)     = percentCorrect(IDX_Test_permuted,Test_labels);

ALGO        = 4;
tic
[~, ClusterCenters]    = kmeans( Train(coreset,:), K );
Timing(ALGO) = toc;
[~,IDX_Train  ]     = pdist2_faster(ClusterCenters,Train,'squaredeuclidean','smallest',1);
[~,IDX_Test   ]     = pdist2_faster(ClusterCenters,Test,'squaredeuclidean','smallest',1);
IDX_Train_permuted  = bestMap( Train_labels, IDX_Train );
IDX_Test_permuted   = bestMap( Test_labels, IDX_Test );
TrainError(ALGO)    = percentCorrect(IDX_Train_permuted,Train_labels);
TestError(ALGO)     = percentCorrect(IDX_Test_permuted,Test_labels);

%% Print out results
for ALGO = 1:4
    fprintf('Training error, %23s: %.2f\n', ALGO_NAMES{ALGO}, TrainError(ALGO) );
end
fprintf('\n');
for ALGO = 1:4
    fprintf('Test error, %23s: %.2f\n', ALGO_NAMES{ALGO}, TestError(ALGO) );
end
fprintf('\n');
for ALGO = 1:2
    fprintf('Timing, %23s: %.2f sec\n', ALGO_NAMES{ALGO}, Timing(ALGO) );
end
for ALGO = 3:4
    fprintf('Timing, %23s: %.2f sec = %.2f + %.2f\n', ALGO_NAMES{ALGO},...
        Timing(2)+Timing(ALGO),Timing(ALGO),Timing(2) );
end
fprintf('Coresets used M=%d (of %d possible, so %.1f%%) points\n', ...
    M, N, 100*M/N );