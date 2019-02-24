%{
 Testing computational speed of various Fast JL transforms
  Feb 11 2019

Y   = S*X, X is M x N and S is m x M

%}
%% Setup paths
% addpath ~/'Google Drive'/GroupDocuments/MatlabUtilities/ % for countSketch
addpath ~/Repos/randomized-algorithm-class/Code/ % from https://github.com/stephenbeckr/randomized-algorithm-class
addpath ~/Repos/hadamard-transform/ % from https://github.com/jeffeverett/hadamard-transform

%% Check Implementations for correctness
X       = randn(2^13,100);
[M,N]   = size(X); m = M/4;

% == Hadamard Code ==
% Check for normalization
Hadamard_teaching_code( eye(4) )
4*fwht(eye(4),[],'hadamard')
%%
tic
Y1 = Hadamard_teaching_code(X);
toc
tic
Y2 = M*fwht(X,[],'hadamard');
toc
tic
Y3 = hadamard_pthreads(X); % my mex code
toc
fprintf('Hadamard code discrepancies: %g and %g\n', norm(Y1-Y2,'fro'), ...
    norm(Y1-Y3,'fro') );

%% == Count sketch --
d       = sign(randn(M,1));
D       = spdiags(d,0,M,M); % bsxfun() is another efficient way to do this
useTranspose    = true;
indx_map        = int64(randi(m,M,1)); % don't do this in C!
Y2 = countSketch_BLAS(X'*D,indx_map,m,useTranspose)';

% Do Count sketch slowly to check
DX      = D*X;
targetRows  = double(indx_map);
Y       = zeros(m,N);
for j = 1:M
    i   = targetRows(j);
    Y(i,:) = Y(i,:) + DX(j,:);
end
fprintf('Count sketch code discrepancies: %g\n', norm( Y - Y2, 'fro' ) )


%% Test speed
N       = 100;
% M_list  = round(logspace( 3, 4, 4 )); % 4 points between 10^3 and 10^4
M_list  = 2.^(10:13);

nTrials     = 10; % get some averages
nAlgos      = 7;
Times       = zeros(nAlgos,length(M_list),nTrials);
% Times_setup = zeros(3,length(M_list),nTrials); % time to make sparse matrix
Times_setup = Times;
ALGO_NAMES = {'Gaussian','FJLT, DCT','FJLT, Hadamard','Count','Very Very sparse',...
    'Very sparse','Sparse'};

for Mi  = 1:length( M_list )
  fprintf('Mi is %d of %d\n', Mi, length(M_list) );
  for trial = 1:nTrials
      
    M   = M_list(Mi);
    m   = round( M/4 );
    
    X       = randn(M,N );
    
    ALGO    = 1; % Gaussian sketch
    tic
    S   = randn(m,M);
    Times_setup(ALGO,Mi,trial)     = toc;
    Y   = S*X; 
    Times(ALGO,Mi,trial) = toc;
    
    ALGO = 2; % Fast JL, DCT
    tic;
    D       = spdiags( sign(randn(M,1)) ,0,M,M);
    ind     = randsample(M,m); % in Stats toolbox
    ind     = randperm(M,m); % faster than randsample, doesn't need toolbox
    Times_setup(ALGO,Mi,trial)     = toc;
    Y       = dct( D*X );
    Y       = Y(ind,:);
    Times(ALGO,Mi,trial) = toc;
    
    ALGO = 3;  % Fast JL, Hadamard
    tic;
    D       = spdiags( sign(randn(M,1)) ,0,M,M);
    %ind     = randsample(M,m); % in Stats toolbox
    ind     = randperm(M,m); % faster than randsample, doesn't need toolbox
    Times_setup(ALGO,Mi,trial)     = toc;
    Y       = hadamard_pthreads( D*X );
    Y       = Y(ind,:);
    Times(ALGO,Mi,trial) = toc;
    
    ALGO = 4; % Count
    tic;
    D               = spdiags( sign(randn(M,1)) ,0,M,M);
    useTranspose    = true;
    indx_map        = int64(randi(m,M,1));
    Times_setup(ALGO,Mi,trial)     = toc;
    Y = countSketch_BLAS(X'*D,indx_map,m,useTranspose)';
    Times(ALGO,Mi,trial) = toc;
    
    
    % Sparse. We can normalize later. Does that help speed?
    for Sparsity = 1:3
        ALGO        = ALGO + 1;
        s           = 2^(1-Sparsity)*sqrt(M);
        density     = 1/(2*s);
        ALGO_NAMES{ALGO} = sprintf('%.1f%% sparse',100*density);
        tic
        S   = sprandn(m,M,density); % this takes longer than the multiply!
        S   = sign(S);
        % SS = logical(S); % alternative
        Times_setup(ALGO,Mi,trial)     = toc;
        Y   = sqrt(s)*(S*X); 
          % is this faster if S is "logical"? that doesn't work,
          % it only has 1 bit, need 2 bits, but sparse of type uint8 not
          % supported
        Times(ALGO,Mi,trial) = toc;
    end
  end
end
%% Plot

Data = Times;
% Data = Times - Times_setup;
mn   = mean(Data,3);

figure(1); clf;
h=loglog( M_list, mn','o-','linewidth',2 );
set(gca,'fontsize',16);
h(2).LineStyle = '--';
h(3).LineStyle = '--';
h(5).LineStyle = ':'; h(6).LineStyle = ':'; h(7).LineStyle = ':';
legend(ALGO_NAMES,'location','northwest','box','off');

% Add something for reference
hold all
y   = M_list/M_list(1);
h2 = loglog( M_list, mn(1)*y.^2, 'k--','DisplayName','M^2','linewidth',2 );
h3 = loglog( M_list, mn(3,1)*...
    (M_list.*log(M_list)/(M_list(1)*log(M_list(1)))),...
    'k-.','DisplayName','M log M','linewidth',2 );
h4 = loglog( M_list, mn(4,1)*y, 'k:','DisplayName','M','linewidth',2 );

xlim([M_list(1),M_list(end)]);
ylim([min(mn(:)),max(mn(:))]);
xlabel('Size M');
ylabel('Time (seconds)');
title('Total times, including setup');
%%
% export_fig 'FastJLtimes_withSetup' '-pdf' -transparent

%%
Data = Times - Times_setup;

mn   = mean(Data,3);

figure(1); clf;
h=loglog( M_list, mn','o-','linewidth',2 );
set(gca,'fontsize',16);
h(2).LineStyle = '--';
h(3).LineStyle = '--';
h(5).LineStyle = ':'; h(6).LineStyle = ':'; h(7).LineStyle = ':';
legend(ALGO_NAMES,'location','northwest');

% Add something for reference
hold all
y   = M_list/M_list(1);
h2 = loglog( M_list, mn(1)*y.^2, 'k--','DisplayName','M^2','linewidth',2 );
h3 = loglog( M_list, mn(3,1)*...
    (M_list.*log(M_list)/(M_list(1)*log(M_list(1)))),...
    'k-.','DisplayName','M log M','linewidth',2 );
h4 = loglog( M_list, mn(4,1)*y, 'k:','DisplayName','M','linewidth',2 );

xlim([M_list(1),M_list(end)]);
ylim([min(mn(:)),max(mn(:))]);
xlabel('Size M');
ylabel('Time (seconds)');
title('Times to apply sketch, excluding one-time setup cost');
%%
% export_fig 'FastJLtimes_excludingSetup' '-pdf' -transparent