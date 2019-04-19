%{
Examples of Localitiy Sensitive Hashing (LSH):

(1) MinHash for Jaccard Similarity
(2) Euclidian Norm
(3) SimHash for Cosine distance

For MinHash, let's suppose we're motivated by hashing documents and then checking for their
 similarity, so we can check if there is plagiarism in a work
See https://www.bowdoin.edu/dean-of-students/judicial-board/academic-honesty-and-plagiarism/examples.html
for plagiarism examples

Stephen Becker
%}

%% Plagiarism Demo: First, convert some sentences/documents to nicer form
sentence = {};
sentence{1} = 'Only two years later, all these friendly Sioux were suddenly plunged into new conditions, including starvation, martial law on all their reservations, and constant urging by their friends and relations to join in warfare against the treacherous government that had kept faith with neither friend nor foe';
sentence{2} = 'Only two years later, all these nice Sioux were suddenly thrust into new types of conditions, including starvation, martial law on all their reservations, and constant urging by their friends and relations to join in warfare against the treacherous government that had kept faith with neither friend nor foe';
sentence{3} = 'In ages which have no record these islands were the home of millions of "Contrast the condition into which all these friendly Indians are suddenly plunged now, with their condition only two years previous: martial law now in force on all their reservations; themselves in danger of starvation, and constantly exposed to the influence of emissaries from their friends and relations, urging them to join in fighting this treacherous government that had kept faith with nobody--neither with friend nor with foe';

% Pre-process, and make k-shingles 
% (and usually you then has the k-shingles down further)
k   = 9;

% https://www.mathworks.com/matlabcentral/answers/3314-hash-function-for-matlab-struct
Engine = java.security.MessageDigest.getInstance('MD5');
for i = 1:3
    s   = sentence{i};
    s   = s( ~isspace( s ) & (s~=',') & (s~='"') ); % remove some characters
    vec = [];
    for start = 1:length(s)-k+1
        ss = s(start:start+k-1); % k-shingle
        % Now, hash it. Use MD5 or SHA since Java does that for me
        % Engine = java.security.MessageDigest.getInstance('MD5');
        Engine = java.security.MessageDigest.getInstance('SHA');
        Engine.update(typecast(uint16(ss), 'uint8'));
        hash    = Engine.digest; % 8 bits per (so 1 byte); keep a few of these
        hash    = uint16(typecast( hash(1:2), 'uint8' )); % and remove signs
%         h       = dec2bin( hash );
%         h       = uint16( bin2dec( h(:)' ) ); % we kept 2 bytes, so 16 bit
        % (Above line is slow, and it merges binary vectors in a funny way.
        % Better is this line below:)
        h       = 2^8*hash(1) + hash(2) + 1; % make it 1-based; be careful to make sure everything is uint16 not uint8 or you have overflow!
        vec     = [vec,h]; % append
        % sprintf('%.2x',double(typecast(hash, 'uint8'))) 
    end
    sentence{i}=vec;
end
% intmax('uint16') % max is 2^6 = 65k
%%
JaccardSim = @(A,B) length(intersect(A,B))/length(union(A,B));
disp( JaccardSim( sentence{1}, sentence{2} ) )
disp( JaccardSim( sentence{1}, sentence{3} ) )
disp( JaccardSim( sentence{2}, sentence{3} ) )
% for i = 1:3, disp(length(sentence{i})); end
%% Apply minhash (many of them), naive version
% Need a universe of all possible entries
% Either take union(...) or use max of uint...
L   = 10; % Number of hashes to draw (L=20 to visualize)
% L   = 1e3; % To check
MinHashSignatures = zeros(3,L);
for ell = 1:L
    P = randperm( intmax('uint16') ); % random permutation of 1, ..., 65k
    for i = 1:3
        MinHashSignatures(i,ell) = min(P(sentence{i}));
    end
end
if L <= 20
    disp(MinHashSignatures)
end
%% Check
for i = 1:3
    for j = (i+1):3
        prob    = sum( MinHashSignatures(i,:) == MinHashSignatures(j,:) )/L;
        fprintf('%d vs %d: JaccDiff is %.2f, %% hash collisions is %.2f\n', ...
            i,j,JaccardSim( sentence{i}, sentence{j} ), prob );
    end
end

%% To tune definition of "neighbors", make bands
% If documents match in *any* band, then declare them
%   a possible neighbor.
rng(0);
b   = 20; % # of bands
r   = 5;  % # hashes per band (if small, then more collisions)
L   = b*r;
MinHashSignatures = zeros(3,b);
for bi = 1:b
    temp = zeros(3,r);
    for ri = 1:r
        P = randperm( intmax('uint16') );
        for i = 1:3
            temp(i,ri) = min(P(sentence{i}));
        end
    end
    % For this band, we have r LSH hashes. Combine these r LSH hashes
    %   by... hashing them together!
    % (This last hash is not a LSH, it's a traditional one)
    for i = 1:3
        Engine = java.security.MessageDigest.getInstance('SHA');
        Engine.update(typecast(uint16(temp(i,:)), 'uint8'));
        hash    = Engine.digest;
        hash    = uint16(typecast( hash(1:2), 'uint8' )); % remove signs
%         h       = dec2bin( hash );
%         h       = uint16( bin2dec( h(:)' ) );
        h       = 2^8*hash(1) + hash(2) + 1; % make it 1-based
        MinHashSignatures(i,bi) = h;
    end
end
MinHashSignatures


%% Try some other hashes, like Euclidean norm distance
% Note: for this LSH, probability of collision isn't identically
% proprtional to the Euclidean distance, but it is a valid LSH
addpath ~/Repos/randomized-algorithm-class/Code/

rng(0);
p   = 100;
N   = 10;
X   = randn(N/2,p);
X   = [ X; X + .1*randn(N/2,p) ]; % so some correlated rows

Dist    = pdist2_faster( X, X,  'sqeuclidean' );

% Now, let's hash these...
a       = .1; 

% Now, combine via banded strategy
rng(0);
b   = 20; % # of bands
r   = 5;  % # hashes per band (if small, then more collisions)
L   = b*r;
Signatures = zeros(N,b);
for bi = 1:b
    temp = zeros(N,r);
    for ri = 1:r
        v   = randn(p,1); v     = v/norm(p); % random unit normal
        bb  = a*rand(1);    % random offset, uniform in [0,a]
        temp(:,ri)  = floor( (X*v + bb)/a );
    end
    % For this band, we have r LSH hashes. Combine these r LSH hashes
    %   by... hashing them together!
    % (This last hash is not a LSH, it's a traditional one)
    for i = 1:N
        Engine = java.security.MessageDigest.getInstance('SHA');
        Engine.update(typecast(uint16(temp(i,:)), 'uint8'));
        hash    = Engine.digest;
        hash    = uint16(typecast( hash(1:2), 'uint8' )); % remove signs
%         h       = dec2bin( hash );
%         h       = uint16( bin2dec( h(:)' ) );
        h       = 2^8*hash(1) + hash(2) + 1; % make it 1-based
        Signatures(i,bi) = h;
    end
end
Signatures


%% Try some other hashes, like SimHash for cosine distances
% For this one, chance of collision is directly proportional to distance

rng(0);
p   = 100;
% N   = 10;
N   = 1e2;
X   = randn(N/2,p);
X   = [ X; X + .1*randn(N/2,p) ]; % so some correlated rows

% Look at cosine distances between all the points in X
nrms    = sqrt( sum(X.^2,2) );
cosDist = real( acos( X*X'./( nrms*nrms' ) ));


% Check if we have collisions at a rate proportional to cosDist: yes!
r   = 1e4;  % repeat it a lot to collect statistics
CollisionFrequency    = zeros(N,N);
temp = zeros(N,1);
for ri = 1:r
    v   = randn(p,1);
    temp  = sign(X*v);
    % This will be slow...
    for i = 1:N
        simInd  = find( temp == temp(i) );
        CollisionFrequency(i,simInd)  = CollisionFrequency(i,simInd) + 1;
    end
end
CollisionFrequency    = CollisionFrequency/r;
TrueFrequency         = 1 - cosDist/pi;
if N <= 10
    disp( CollisionFrequency )
    disp( TrueFrequency )
    disp( CollisionFrequency - TrueFrequency )
else
    [TrueFreq_sorted, sort_ind] = sort( TrueFrequency(:) );
    figure(1); clf;
    scatter( TrueFreq_sorted(:), CollisionFrequency( sort_ind ),'r.' );
    hold all
    line( [0,1],[0,1],'linestyle','--','color','k')
    xlabel('True cosine distance');
    ylabel('Frequency of LSH collision');
end

% We can also combine these in the same banding technique...