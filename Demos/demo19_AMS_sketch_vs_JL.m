%{
Short demonstration of the AMS sketch
AMS is named for the authors of this paper:
 "The space complexity of approximating the frequency moments" (N. Alon, Y. Matias, and M. Szegedy, STOC 1996)

There have been many updates, improvements to this sketch over the years.
I'm following *some* of the improvements, mainly as described
in Cormode's 2013 lecture http://dimacs.rutgers.edu/~graham/pubs/html/TalkSimons13.html
but see also Cormode's 2011 monograph for more precise statements,
 "Sketch techniques for approximate query processing" 
  Foundations and Trends in Databases
  http://www.cs.umass.edu/~mcgregor/711S12/sketches1.pd


Stephen Becker
This code uses the AMS_sketch.m code, in the ../Code subdirectory
    on this same github repo.
%}
rng(0);
addpath ~/Repos/randomized-algorithm-class/Code/

p   = 1e4;
n   = 1e2;
w   = 2^10;
d   = 7;

X   = randn(p,n);
% Y   = randn(p,n);
% C   = AMS_sketch( X, w, d );    % calls several count sketches
% CY  = AMS_sketch( Y, w, d );
% C2  = AMS_sketch( X - 3.14*Y, w, d );
% norm( (C-3.14*CY) - C2, 'fro' ) % confirm linearity of sketch

if false
    addpath('~/Google Drive/TeachingDocs/APPM4720_5720_Spring19_Randomized/Code');
    load mnist_data_all
    X   = Train';
    Xt  = Train;
    d   = 7;
    m   = 20;
    n   = size(X,2);
end

columnNorms     = @(X) sqrt( sum(X.^2) );
%%
saltSeed    = 1;
tic
if n > 1e4
    C   = AMS_sketch( Xt, w, d, 'transposedX', true ); % calls several count sketches
else
    C   = AMS_sketch( X, w, d ); % calls several count sketches
end
toc
%%
cNorms  = zeros(d,n);
for k = 1:d
    CC  = C( (1+(k-1)*w):k*w, :);
    cNorms(k,:)     = columnNorms( CC );
end
cNormEstimate   = median( cNorms, 1 );
cNormEstimate_variant   = sqrt( mean( cNorms.^2, 1 ) ); % e.g., w<-- w*d, d<-- 1
%%
figure(1); clf;
scatter( columnNorms(X), cNormEstimate, 'o' )
hold all
scatter( columnNorms(X), cNormEstimate_variant, 'x' )
line( [96,104],[96,104] );
axis equal % make it dramatic!
legend('Median of rows','Mean of rows');
%%
figure(1); clf;
histogram( cNormEstimate./columnNorms(X), 'binwidth',.01 )
hold all
histogram( cNormEstimate_variant./columnNorms(X), 'binwidth',.01 )
legend('Median of rows','Mean of rows');

