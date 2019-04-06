% demonstrate Compressed Sensing ideas
% This code requires CVX (cvxr.org)

rng(0);

N       = 100; % dimensionality of signal
s       = 5;   % sparsity of signal

x0      = zeros(N,1);
x0( randperm(N,s) )     = rand(s,1); % random entries

% Try this for different values of m. How low can you go?
% m       = 4*s;
m       = round( 2.5*s ); % theoretical lower limit is 2*s
A       = randn(m,N);   % Sensing matrix

% figure(1); clf; imagesc(A); axis image

y       = A*x0;

cvx_begin quiet
    variable x1(N)
    minimize norm(x1,1)
    subject to
        A*x1 == y
cvx_end

x1( abs(x1) < 1e-9 ) = 0;

x2  = pinv(A)*y; % least-squares solution

%% Plot
figure(1); clf;
stem( find(x0), x0(find(x0)), 'd' , 'markersize',10);
hold all
stem( find(x1), x1(find(x1)), 'o','MarkerFaceColor','r');
stem( x2, '*' );
set(gca,'fontsize',16)
legend('Original','l1','l2','location','best');

%% Can we do this with a combinatorial algorithm? No
% Assuming we know s
nchoosek(N,s)  % # of permutations to try

% Make a list of all permutations... or not. Pretty slow!
tic;
list = nchoosek( 1:N, 2 );
toc
tic;
list = nchoosek( 1:N, 3 );
toc
tic;
list = nchoosek( 1:N, 4 );
toc