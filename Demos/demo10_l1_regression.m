%{

(1) Motivate l1 or l_p (1<=p<=infty) regression

(2) Sketching approaches
Ref: "The Fast Cauchy Transform..." by Clarkson et all 2016
in SIAM J Sci Comp, ?http://epubs.siam.org/doi/10.1137/140963698

see also David Woodruff's 2014 monograph

Stephen Becker

%}

%% Why l1 regression? Let's do a 1D example
rng(0);
x   = 2;
A   = (1:6)';
[M,N]   = size(A); % N = 1 since 1D example
z   = .3*randn(M,1);
b   = A*x + z;

b(6) = 1.2; % gross corruption, e.g., someone typed "1.2" instead of "12"

% Notation is funny: "x" is really like a slope,
%   and "A" is really like "x"
figure(1); clf;
h=plot( A, b, 'd', 'DisplayName', 'Samples','markersize',10,...
    'markerfacecolor','b');
hold all
% plot( A, A*x, '-', 'DisplayName', 'True data','linewidth',2);

% Find line of best fit, in l2 sense
xLS     = A\b;
plot( A, A*xLS, '--', 'DisplayName', 'l2 fit','linewidth',2);
cvx_begin quiet
    variable x1(N)
    minimize norm( A*x1 - b , 1 )
cvx_end
plot( A, A*x1, '--', 'DisplayName', 'l1 fit','linewidth',2);
legend('location','northwest');
set(gca,'fontsize',18);


%% Check basic Johnson-Lindenstrauss results: preserve distance in lp sense
rng(0);
% nPoints     = 1e3;
M           = 5e2; % dimension of each point
% make the data points, some of them sparse, some weird distributions, ...
A   = [log(abs(gallery('binomial',M))), gallery('dramadah',M), ...
    gallery('cauchy',M), gallery('hanowa',M), gallery('lotkin',M) ];
nPoints     = size(A,2);

clf; cspy(A); title('Depiction of matrix "A"');
% normalize it
nrms2_A = sqrt( sum(A.^2,1) );
nrms1_A = sum(abs(A),1);
% A       = bsxfun( @times, A, 1./nrms );

%% Take the l2 sketch
m   = round(.3*M);
addpath ~/Repos/randomized-algorithm-class/Code/
%% Gaussian sketch
S   = randn(m,M)/sqrt(m);
SA  = S*A;
%% FJLT sketch
ind = randperm(M,m);
SA  = dct(spdiags(sign(randn(M,1)),0,M,M)*A);
SA  = sqrt(M/m)*SA(ind,:);
%% Cauchy sketch
% Same as student-t with 1 degree of freedom
S   = trnd(1,m,M);
SA  = 1/sqrt(m)*S*A;
%% Check if we've preserved l1 and l2 norms
nrms2    = sqrt( sum(SA.^2,1) );
nrms1    = sum(  abs(SA),1);
figure(1); clf;
subplot(1,2,1);
histogram( nrms2./nrms2_A,'Normalization','probability' )
xlim([0,2]);
title('$\|Sx\|_2/\|x\|_2$','interpreter','latex','fontsize',20);

subplot(1,2,2);
histogram( nrms1./nrms1_A,'Normalization','probability' )
title('$\|Sx\|_1/\|x\|_1$','interpreter','latex','fontsize',20);

%% Zoom in on histograms for the case of Cauchy sketch
figure(1); clf;
BMIN    = .5;
BMAX    = 40;
subplot(1,2,1);
histogram( nrms2./nrms2_A,'BinLimits',[BMIN,BMAX] ,'Normalization','pdf')
% xlim([0,2]);
title('$\|Sx\|_2/\|x\|_2$','interpreter','latex','fontsize',20);

subplot(1,2,2);
histogram( nrms1./nrms1_A ,'BinLimits',[BMIN,BMAX],'Normalization','pdf')
title('$\|Sx\|_1/\|x\|_1$','interpreter','latex','fontsize',20);



%% Interpret another way: use p-stable to estimate p-norms
innerProds    = sqrt( abs(sum(SA,1)) ); % abs after sum

figure(1); clf;
BMIN    = .5;
BMAX    = 10;
histogram( innerProds./nrms1_A,'BinLimits',[BMIN,BMAX] ,'Normalization','pdf')
% xlim([0,2]);
title('$\sqrt{E|\langle x, s \rangle|^2} /\|x\|_1$','interpreter','latex','fontsize',20);






%% Regression
rng(0);
M   = 1e3;
N   = 1e2;
A   = rand(M,N);
x0  = randn(N,1);
b   = A*x0 + randn(M,1);

% Solve large problem for reference solution
tic
cvx_begin quiet
  variable x(N)
  minimize norm( A*x - b, 1 )
cvx_end
toc
xRef    = x;
%% Make well-conditioned basis
rng(1);
m   = round(.5*M);
for i = 1:2
    switch i
        case 1
            % Cauchy sketch:
            S   = trnd(1,m,M)/sqrt(m);
            fprintf('\nUsing Cauchy sketch\n');
        case 2
            % Gaussian sketch
            S   = randn(m,M)/sqrt(m);
            fprintf('\nUsing Gaussian sketch\n');
    end

SA  = S*A;


[Q,R]   = qr(SA,0);
Q   = A/R;
% estimate l1 leverage scores
levScores   = sum( abs(Q), 2 );
% weighted sampling
ind         = randsample(M,round(.5*M), true,levScores );

%  == Solve smaller problem
tic
cvx_begin quiet
  variable x(N)
  minimize norm( A(ind,:)*x - b(ind), 1 )
cvx_end
toc
er1=norm( x - xRef )/norm(xRef );
er2=norm( A*x - b, 1 )/norm( A*xRef - b, 1 ) - 1;
fprintf('||x-xRef| is %.2e, ||Ax-b||_1/||AxRef-b||_1-1 is %.2e\n', er1,er2);

end