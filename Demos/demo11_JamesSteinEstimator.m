%{

James-Stein Estimate
  Proof that the MLE is not admissible (in dimensions p>=3 at least)

see
http://statweb.stanford.edu/~ckirby/brad/LSI/chapter1.pdf
and
http://statweb.stanford.edu/~ckirby/brad/other/CASI_Chap7_Nov2014.pdf
for connection to Empirical Bayes...

 (it's also similar to the idea of control variates)

%}

p       = 50; % dimension
sigma   = .3;

nReps   = 1e3;

mu      = zeros(p,1);
v       = .5*ones(p,1); % arbitrary fixed vector
% v       = .1*randn(p,1);

sampleMeans     = zeros(p,2);
firstCoordinate = zeros(nReps,2);
avgError        = zeros(nReps,2);
errFcn          = @(xhat) norm(xhat-mu)^2;
for r = 1:nReps
    
    y   = mu + sigma*randn(p,1);
    
    % MLE is y
    sampleMeans(:,1)    = sampleMeans(:,1) + y;
    firstCoordinate(r,1)= y(1);
    avgError(r,1)       = errFcn(y);
    
    % James-Stein estimator
    xhat    = (1 - (p-3)*sigma^2/( norm(y-v)^2 ) )*(y-v) + v;
    sampleMeans(:,2)    = sampleMeans(:,2) + xhat;
    firstCoordinate(r,2)=xhat(1);
    avgError(r,2)       = errFcn(xhat);
    
end
sampleMeans = sampleMeans/nReps;

%% Analyze results
figure(1); clf;
boxplot( avgError,'Labels',{'MLE','James-Stein'} )
set(gca,'fontsize',18);
title('Values of $\|\hat{x} - \mu\|_2^2$','interpreter','latex')

%% Look at estimate of the first coordinate: is it biased?
% (mu and v are all the same in all coordinates, so just pick the first
%  coordinate, since then it's easy to show graphically)
figure(1); clf;
boxplot( firstCoordinate,'Labels',{'MLE (unbiased)','James-Stein (biased!)'} )
set(gca,'fontsize',18);
title('First coordinate of the estimate');
line([-.5,2.5],[0,0],'color','k','linestyle','--')