%{
Demo for 1D leverages scores

Fig. 2 in Section 6.1 of 
"?Randomized algorithms for matrices and data"
(Mahoney, 2011, ?http://arxiv.org/abs/1104.5557 )

is misleading, since it discussing perturbing regressors and data

Consider a 1D regression problem,
min_{beta} || x*beta - y ||_2

where x is a n x 1 vector of covariates, and y are the data.

By perturbing x, we change the solution.
This is captured by the leverages scores of x.
Since x is a vector, leverage scores are simply proportional
to the magnitude of each entry.

Therefore, the idea of "leverage" is that if we perturb
entries of x that have more leverage, i.e., that are large
in magnitude, then the effect on the regression is greater.

Do we observate that?

Stephen Becker, Feb 14 2018

%}

rng(0);
n   = 5;
% x   = 1:n; 
% x   = -n:-1;
x   = -2:2;


x = x';
slope   = 1;
y   = slope*x + .1*randn(n,1);

slopeEst    = x\y;

delta   = -.9; % try +/-
i = 1;
figure(1); clf;

for i = 1:5
    
    subplot(2,3,i);
    
    plot( x, y ,'o','markersize',10)
    hold all
    plot( x, slopeEst*x, '--' ,'linewidth',2)
    xx = x; xx(i) = xx(i) + delta;
    plot( xx(i), y(i),'s','markersize',10,'MarkerFaceColor','k');
    line( [x(i),xx(i)], y(i)*[1,1]);
    slopeEst_perturbed    = xx\y;
    plot( x, slopeEst_perturbed*x, '-','linewidth',2 )
    xlim([min(x)-abs(delta),max(x)+abs(delta)]);
    title(sprintf('Moving %d^{th} data point',i))

end
