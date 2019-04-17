%{
unweighted k-Nearest Neighbors

Note: everything would probably be faster if we stored the data
 where columns are new data points, not rows, but we're sticking
 with the row convention since it's more common (and Matlab uses it
 for their functions, even though it goes against their natural
 data structure).

Compare also with Matlab's knnsearch implementation
  If dimension p < 10 then this can exploit a kd-tree at training time,
  but the complexity of that scales very poorly with dimension,
  so not applicable to MNIST without doing some dimensionality reduction.

Stephen Becker
%}
clc

addpath ~/Repos/randomized-algorithm-class/Code/
% Load MNIST data:
addpath('~/Google Drive/TeachingDocs/APPM4720_5720_Spring19_Randomized/Code');
load mnist_data_all

percentCorrect = @(labels1,labels2) length(find(labels1(:)==labels2(:)))/length(labels1(:));

rng(1);
% Try with various sizes to get an idea how it scales...
% test_subset     = randsample( 1e4, 1e3  );
test_subset     = randsample( 1e4, 1e2  );

TestPoints      = Test(test_subset,:);
TestLabels      = Test_labels( test_subset );
K               = 10;   % # of nearest neighbors to use
%% Try k-NN using true distances, so plain implementation

fprintf('\n-- Vanilla k-NN\n  Finding pairwise distances\n');
t1  = tic;
tic
% ind   = dsearchn( Train, TestPoints ); % Slow
% D     = pdist2( Train, TestPoints ); % Slow
% for 1e3 test points,takes 34 sec with pdist, or 2.87 with pdist_faster
D   = pdist2_faster( Train, TestPoints );
toc
fprintf('  Sorting those distances\n');
tic;
[~,Ind]   = sort(D); % per row, sort the columns
toc

fprintf('  Final processing\n'); % Find the labels of the neigbhors
labels      = Train_labels( Ind(1:K,:) );
prediction  = mode( labels, 1 );
pc          = percentCorrect( prediction, TestLabels );
fprintf('  Standard k-NN has %.1f%% accuracy\n', pc*100 );

t_plain     = toc(t1);



%% Do with Matlab's knn
% If p >= 10, it won't se a kd tree (see KDTreeSearcher)
% See "Classificiation Using Nearest Neighbors" help topic
% Mdl     = ExhaustiveSearcher( Train, 'Distance', 'seuclidean' );
% [ind2,dist_ind] = knnsearch(Mdl,TestPoints,'k',K); % gave NaNs
if size(TestPoints,1) < 500
    fprintf('\n-- Vanilla k-NN via Matlab''s implementation\n');
    t1 = tic;
    [idx,d]     = knnsearch( Train, TestPoints, 'K', K );
    toc(t1)
    
    labels  = zeros(size(TestLabels));
    t2=tic;
    for i = 1:size(TestPoints,1)
        labels(i)  = mode(Train_labels( idx(i,:) ));
    end
    toc(t2)
    t_Matlab = toc(t1);
    pc_Matlab          = percentCorrect( labels, TestLabels );
    fprintf('  k-NN via Matlab has %.1f%% accuracy\n', pc_Matlab*100 );
else
    pc_Matlab = nan;
    t_Matlab = Inf;
end


%% Do it with LSH and bands
% b = 1, r = 30 is very bad (very few neighbors, but lots of false
% negatives); better to increase b and decrease r
% (for Cosine vs Euclidean distances, parameters will vary
%  and for Euclidean distance, also work with "a" parameter)
b       = 15; % number of bands (decrease this to reduce # neighbors found)
r       = 3;  % hashes per band (increase this to reduce # neighbors found)
a       = 5e2; % controls fineness; check length( unique( Train_hashed(:,1) ) )

t1  = tic;
rng(1);
p       = size(Test,2);
fprintf('\n-- k-NN via LSH\n');
tic

COSINE_DISTANCE = false;
neighborList    = zeros( size(Train,1), length(TestLabels),'logical' );

for bi = 1:b
    if COSINE_DISTANCE
        % Cosine distance "SimHash"
        Omega   = randn(p,r);
        LSH     = @(X) sign((X-mean(X,2))*Omega);
        
        Train_hashed     = LSH(Train);
        Test_hashed      = LSH(TestPoints);
        
        innProd       = Train_hashed*Test_hashed';
        neighborList  = neighborList | (innProd==r); % binary "or"
    else
        % Euclidean distance hash
        V   = randn(p,r);
        bb  = a*rand(1,r);
        LSH = @(X) floor( (X*V + bb )/a );
        Train_hashed     = LSH(Train);
        Test_hashed      = LSH(TestPoints);

        % Do group updates, assuming we only have a few hash values
        universe = unique( Test_hashed );
        tempList = zeros( size( neighborList), 'uint8' );
        % Avoid "unique" call by using big matrix
        for ri = 1:r
            for val_i = 1:length(universe)
                val     = universe( val_i ); % bucket value
                ind_test    = find( Test_hashed(:,ri) == val );
                ind_train   = find( Train_hashed(:,ri) == val );
                tempList( ind_train, ind_test ) = tempList(ind_train, ind_test) + 1;
            end
        end
        % Need to hash all ri things together (need them *all* to agree)
        %  or, check when tempList == r
        neighborList = neighborList | (tempList==r);
        
    end
end
toc

fprintf('  Reduced # of neighbors to %.1f%%\n', 100*nnz(neighborList)/numel(neighborList) );

labels  = zeros(size(TestLabels));
tic
for i = 1:size(TestPoints,1)
    ind     = find( neighborList(:,i) );
    [~,ind2] = sort( pdist2_faster( Train(ind,:), TestPoints(i,:) ) );
    KK       = min( length(ind2),K );
    labels(i)  = mode(Train_labels( ind(ind2(1:KK)) ));
end
toc
pc_LSH          = percentCorrect( labels, TestLabels );
fprintf('  LSH k-NN has %.1f%% accuracy\n', pc_LSH*100 );

t_LSH   = toc(t1);
%% Overall
fprintf('\n== SUMMARY ==\n %6d training points, %5.1f s via plain k-NN (%5.1f via Matlab''s), %5.1f s via LSH k-NN\n',...
    size(TestPoints,1), t_plain,t_Matlab, t_LSH );
fprintf('\tand respective accuracies: %.1f%%, %.1f%% and %.1f%%\n', pc*100, pc_Matlab*100, pc_LSH*100 );

