%{
With what frequency do names occur in the US?
Use social security data

From https://www.ssa.gov/oact/babynames/limits.html
download https://www.ssa.gov/oact/babynames/names.zip
(about 9 MB)

Top 10 from 2017 are here: https://www.ssa.gov/oact/babynames/


This demo applies the CountMin sketch to estimate the frequency
 of occurrence of each name (using less memory than the straightforward
 data structure).
 Note: for both the "straightforward" data structure (which we call
 "fullData") as well as the sketch, we're really storing a hash table (I
 use the SHA hash, since it's easy to call via Matlab/Java),
 so I don't actually store the names themselves! But if you think of a
 name, we can then hash it, and check if it's in the table.

 Note: for the hash table, I use the first 2 bytes of the SHA hash,
 so about 65k unique buckets. There are at least 50,686 names in the
 database, so about 77% of buckets are occupied, so there are a lot
 of name collisions. This is bad! The fix is to use more bytes
 of the SHA hash, but I'm too lazy to implement that now.
 [Update: I fixed my laziness and used 2^20 buckets]

 Some of the hashes, the ones used in the sketch, only need pairwise
 independence, and you can do things faster than SHA and MD5
 e.g., if w is a prime number, then drawing a and b randomly
    from [0, w-1], the function h(x) = a*x + b (mod w)
    has the pairwise independence probability: the chance
    of two inputs colliding (over the randomness of choosing a and b)
    is 1/w. 

Stephen Becker, University of Colorado

Reference: 
Graham Cormode, ?http://dimacs.rutgers.edu/~graham/pubs/html/TalkSimons13.html 
and his monograph: "?Sketch techniques for approximate query processing"
2011 (?http://www.cs.umass.edu/~mcgregor/711S12/sketches1.pdf)
%}

% nNameBuckets = intmax('uint16'); % too small
% nNameBuckets = intmax('uint32'); % too big
nNameBuckets    = 2^20;
fullData     = zeros(nNameBuckets,1);
fullDataNames     = cell(nNameBuckets,1);
if 2==exist('demo18_data.mat','file') %&& false
    load demo18_data
    fullData = full( fullDataSparse );
else
    tic
    fprintf('Reading in year     ');
    for yr = 1880:2017
        fprintf('\b\b\b\b%d',yr);
        prfx    = '~/Downloads/names';
        filename = fullfile(prfx, sprintf('yob%d.txt',yr) );
        fid     = fopen(filename);
        data    = textscan( fid, '%s%c%d','Delimiter',',');
        names   = data{1}; % data{2} is gender, 'M' or 'F'
        occurences  = data{3};
        for line = 1:length(names)
            Engine  = java.security.MessageDigest.getInstance('SHA');
            name    = lower( names{line} );
            Engine.update(typecast(uint16(name), 'uint8'));
            hash    = Engine.digest; % 8 bits per (so 1 byte); keep a few of these
            hash    = uint32(typecast( hash(1:3), 'uint8' )); % and remove signs
            %         h       = dec2bin( hash );
            %         % hash is 0 to nNameBuckets-1, so need a +1 offset
            %         h       = uint16( bin2dec( h(:)' ) ) + 1;  % NO, not quite
            %         right...
            h       = 2^16*hash(1) + 2^8*hash(2) + hash(3);
            h       = mod( h, nNameBuckets ) + 1;
            
            fullData(h) = fullData(h) + occurences(line);
            
            % Also, add the name to the list of names
            if isempty( fullDataNames{h} )
                fullDataNames{h} = name;
            elseif isempty( strfind(fullDataNames{h},name) )
                fullDataNames{h} = [fullDataNames{h},',',name];
            end
                
            
        end
        
        fclose(fid);
    end
    fprintf('  finished.\n'); % take about 3.8 minutes
    toc
    fullDataSparse = sparse( fullData ); % compress from 8 MB to 1.4 MB; fullDataNames is 18 MB
    save demo18_data fullDataSparse fullDataNames % .mat file compresses it anyhow...
end
%%
fprintf('Found at least %d distinct names (maybe more, since could be collisions)\n', nnz(fullData) );
fprintf(' And there were %.1f million people in the dataset\n', sum(fullData)/1e6 );
% Find collisions this way:
fullDataCollisions     = zeros(nNameBuckets,1);
for j = 1:nNameBuckets
    if ~isempty( fullDataNames{j} )
        str = fullDataNames{j};
        fullDataCollisions(j) = length( strfind(str,',') ) + 1;
    end
end
numberUniqueNames = sum( fullDataCollisions );
fprintf(' Checking for collisions, we found exactly %d distinct names, so %d collisions\n', numberUniqueNames, numberUniqueNames-nnz(fullData) );
collisionIndex = find( fullDataCollisions > 1 );
fprintf(' For example, a few collisions:\n');
fullDataNames{ collisionIndex(1:2) }

%% Warning...
% Lot's of bugs because of datatypes, e.g.,
% j = 981698; 
% typecast(j, 'uint8')
% typecast( uint32(j), 'uint8' )
% The above two things are NOT the same!!

% Most bugs are hopefully fixed!

%% Try CountMin sketch
% We could have applied this as we read in the data files, since
%   it is a linear sketch, so easy to update. But since we have
%   the full data anyhow, let's do it the easy way.
d   = 7; % e.g., ceil( log2( 1/.01 ) ), so result holds with 99% chance;
w   = 2^8;  % number of buckets per each hash
C   = zeros( d, w );
 
for j = find( fullData )' % only loop over non-empty ones rather than for j = 1:nNameBuckets
    Engine  = java.security.MessageDigest.getInstance('SHA');
    Engine.update(typecast(uint32(j), 'uint8')); % uint16(j) is a bug (if j is too big)
    L       = typecast( Engine.digest, 'uint8' );  % make it non-negative
    for k = 1:d
        ell     = L(k) + 1; % make it 1-based not 0-based indexing
        C( k, ell ) = C( k, ell ) + fullData( j );
    end
end

%% And repeat, but with more buckets
d2   = 7; 
w2  = 2^12;
C2  = zeros( d2, w2 );
for j = find( fullData )'
    Engine  = java.security.MessageDigest.getInstance('SHA');
    Engine.update(typecast(uint32(j), 'uint8'));
    L       = uint16(typecast( Engine.digest, 'uint8' ));
    for k = 1:d2
        ell     = L(2*k-1)*2^8 + L(2*k);
        ell     = mod( ell, w2 ) + 1;
        C2( k, ell ) = C2( k, ell ) + fullData( j );
    end
end
%% Difference in sizes:
fprintf('Full data has %d entries\n', length(fullData) );
fprintf('  CountMin structure has %d = %d x %d entries, so %.1fx compression\n', ...
    d*w, d, w, length(fullData)/(d*w) );
fprintf('  and more accurate CountMin structure has %d = %d x %d entries, so %.1fx compression\n', ...
    d*w2, d2, w2, length(fullData)/(d2*w2) );

% More accurate estimate of compression ratio
fullDataSparse = sparse( fullData );
stat    = whos('fullData'); b1  = stat.bytes;
stat    = whos('fullDataSparse'); b2  = stat.bytes;
stat    = whos('C'); b3  = stat.bytes;
stat    = whos('C2'); b4  = stat.bytes;
kB      = 1/1024;
fprintf('Naive: %.1f kB, compressed naive: %.1f kB, CountMin: %.1f kB, CountMin v2: %.1f kB\n', ...
    b1*kB, b2*kB, b3*kB, b4*kB );
fprintf('  So compression ratios %.1fx, %.1fx, %.1fx, %.1fx (relative to naive)\n', ...
    b1/b1, b1/b2, b1/b3, b1/b4 );
fprintf('  ...compression ratios %.1fx, %.1fx, %.1fx, %.1fx (relative to compressed naive)\n', ...
    b2/b1, b2/b2, b2/b3, b2/b4 );
%% Now, try it out
name    = 'james'; % Most popular name in database
% name    = 'alexander';
% name    = 'samantha';
% name    = 'john';
% name    = 'hendrix'; % works
% name    = 'sophie';
% name    = 'sophia'; 
% name    = 'marta';
% name    = 'abigayll';
% name    = 'Isabella';
% name    = 'ryan';
% name    = 'padraig';
% name    = 'kathryn'; % known collision
% name    = 'fortino';
% name    = 'tomasz'; % known collision
% name    = 'stephen';

totalNames  = sum( fullData );
%totalNames  = sum( C(:) )/d; % equivalent
% sum(C,2) - totalNames  % sanity check for debugging purposes

% Figure out the index j
Engine  = java.security.MessageDigest.getInstance('SHA');
Engine.update(typecast(uint16(lower(name)), 'uint8'));
hash    = Engine.digest; % 8 bits per (so 1 byte); keep a few of these
hash    = uint32(typecast( hash(1:3), 'uint8' )); % and remove signs
j       = mod( 2^16*hash(1) + 2^8*hash(2) + hash(3), nNameBuckets ) + 1;



% And now try the CountMin sketch
Engine  = java.security.MessageDigest.getInstance('SHA');
Engine.update(typecast(uint32(j), 'uint8'));
L       = uint16(typecast( Engine.digest, 'uint8' ));
c       = Inf;
c2      = Inf;
for k = 1:d
    ell     = L(k) + 1; % make it 1-based not 0-based indexing
    c       = min( [c, C( k, ell )] );
end
for k = 1:d2
    ell2    = L(2*k-1)*2^8 + L(2*k);
    ell2    = mod( ell2, w2 ) + 1;
    c2      = min( [c2, C2( k, ell2 )] );
end
fprintf('\nName: %s\n', name);
fprintf('True frequency is\t%.5f%%\n', 100*fullData(j)/totalNames );
fprintf('Estimated frequency is\t%.5f%% (with CountMin sketch)\n', 100*c/totalNames );
fprintf('Estimated frequency is\t%.5f%% (with larger CountMin sketch)\n', 100*c2/totalNames );
if fullDataCollisions(j) > 1
    fprintf(' Careful! There were other names with hash collisions:  ');
    disp( fullDataNames{j} );
elseif fullDataCollisions(j) == 0
    fprintf(' Careful! This name was not in database, we''re getting only noise\n');
end