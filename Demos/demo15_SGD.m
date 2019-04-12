%{
Stochastic Gradient Descent (SGD) and variants,
 demonstrated on the primal support vector machine (SVM) problem
 with MNIST data.

Primal SVM:

min_w  ||w||^2/2 + C \sum_i^N hinge( y_i w'*x_i )
where
hinge(a) = max( 0, 1 - a )

We will pick 2 classes from MNIST (not all 10)
If we pick class 0 and class 1, it's very easy (accuracy is > 99% in 1
step)
Try harder classes, like 4 and 9

Compare SGD with batch gradient descent,
 as well as
 - minibatch
 - SAGA (minibatch)
 - SVRG (minibatch)

It's all quite messy, since stepsizes are a big issue (and other
parameters)

%}

addpath ~/Repos/randomized-algorithm-class/Code/
addpath('~/Google Drive/TeachingDocs/APPM4720_5720_Spring19_Randomized/Code');

load mnist_data_all
percentCorrect = @(labels1,labels2) length(find(labels1==labels2))/length(labels1);

%% In prep for SVM, reduce to just two labels
A   = 4; B   = 9;       % harder to distinguish
% A   = 0; B   = 1;     % very easy to distinguish

indx   = (Train_labels==A) | (Train_labels==B);
% convert to -1, +1 labels
Train_labels_2class     = Train_labels( indx );
Train_labels_2class( Train_labels_2class==A ) = -1;
Train_labels_2class( Train_labels_2class==B ) = +1;
Train_2class            = Train( indx, : );

indx   = (Test_labels==A) | (Test_labels==B);
Test_labels_2class     = Test_labels( indx );
Test_labels_2class( Test_labels_2class==A ) = -1;
Test_labels_2class( Test_labels_2class==B ) = +1;
Test_2class            = Test( indx, : );

hist(Train_labels_2class ) % make sure it looks OK, about equal
%% Plot hinge loss
x   = linspace(-3,3,40);
plot( x, max(0,1-x),'linewidth',2)
line([-3,3],[0,0],'linestyle','--','color','k');
line([0,0],[-1,4],'linestyle','--','color','k');
ylim([-1,4]); title('Hinge loss');

%% Try SVM... but apply to just two classes
% Pick the dataset
X       = Train_2class;
X       = [X,ones(size(X,1),1)]; % allow for an offset
y       = Train_labels_2class; % this is now -1, +1

yX      = bsxfun( @times, X, y );
[N,p]   = size( X );

C1      = 1e-2;  % constant for SVM model
C2      = 1/N;   % constant for SVM model

maxIts  = 5e2;
minibatch_ratio     = .1; % 10% sampling
minibatch_n         = round(minibatch_ratio*N);

ALGONAMES = {'Gradient Descent','SGD fixed stepsize','SGD decaying stepsize',...
    'SGD minibatch','SAGA minibatch','SVRG'};
errList       = zeros(length(ALGONAMES),maxIts,3); % 3 types of errors
for ALGO = 1:length(ALGONAMES)
    fprintf('\nAlgorithm: %s\n', ALGONAMES{ALGO} );
    
    w       = zeros(p,1); % our variable
    switch ALGO
        case {1,2}
            decay_gamma     = false;
        otherwise
            decay_gamma     = true;
    end
    switch ALGO % choose learning rate
        case {1,4}
            gamma   = 1e-5; % 1e-2 is too big
        case 5 % SAGA
            gamma   = 1e-4; % do it with minibatch
        case 6 % SVRG
            gamma   = 1e-6;
        otherwise
            gamma   = 1e-5;  % SGD needs smaller stepsize
    end
    
    for k = 1:maxIts
    
        %  f(w) = C1 ||w||^2/2 + C2 ones(n,1)*hinge( diag(y)*X*w )
        % where
        % hinge(a) = max( 0, 1 - a )
        % so d(hinge)/da = { -1 (a <= 1); 0 (a > 1) }

        switch ALGO
            case 1
                a   = yX*w;     % helper variable
                grad  = yX'*( -(a<=1) ); % full gradient step
            case {2,3}
                ind     = randperm(N,1);
                a       = yX(ind,:)*w;     % helper variable
                grad    = N*yX(ind,:)'*( -(a<=1) );
            case {4} % minibatch
                ind     = randperm(N,minibatch_n);
                a       = yX(ind,:)*w;     % helper variable
                grad      = (N/minibatch_n)*yX(ind,:)'*( -(a<=1) );
            case {5} % SAGA
                if k==1
                    % First iteration is special: make full pass through
                    % data
                    a       = yX*w;     % helper variable
                    grad    = yX'*( -(a<=1) ); % full gradient step
                    a_storage   = a; % store this
                    grad_storage= grad;
                else
                    %ind     = randperm(N,1);
                    ind     = randperm(N,minibatch_n);
                    a       = yX(ind,:)*w;     % helper variable
                    grad_ind_new    = yX(ind,:)'*( -(a<=1) );
                    
                    % Combine:
                    grad_ind_old  = yX(ind,:)'*( -(a_storage(ind)<=1) );
                    grad    = N/minibatch_n*grad_ind_new - ...
                        N/minibatch_n*grad_ind_old + grad_storage;
                    % Update storage table:
                    a_storage(ind)      = a;
                    grad_storage        = grad_storage ...
                        - grad_ind_old + grad_ind_new;
                end
            case 6
                % SRVG                
                a   = yX*w;     % helper variable
                grad  = yX'*( -(a<=1) ); % full gradient step
                % Make a bunch of micro steps now
                z   = w;
                for kk = 1:50
                    ind = randperm(N,round(N/minibatch_n));
                    a_z     = yX(ind,:)*z;     % helper variable
                    grad_z    = yX(ind,:)'*( -(a_z<=1) );
                    z       = z - gamma*( C1*z + 1/minibatch_n*grad_z +C2*grad );
                end
                w   = z;
                
        end
        
        if decay_gamma && ~mod( k, 50 )
            gamma   = gamma/2;
        end
    
        % Combine to get full gradient, take gradient descent step
        if ALGO ~= 6 % SVRG does its own update
            w   = w - gamma*(C1*w + C2*grad );
        end
        
        % Record metrics:
        % Cost function (expensive to calculate... for academic purposes)
        Xw  = X*w;
        f   = C1*norm(w)^2/2 + C2*sum( max(0,1-y.*Xw) );
        % Percent correct (pc) for test/train
        IDX_Train   = sign( Xw ); % no need to find best permutation
        pc          = percentCorrect(IDX_Train,y);
        IDX_Test    = sign( Test_2class*w(1:end-1) + w(end) ); % allow offset
        pc_test     = percentCorrect( IDX_Test, Test_labels_2class); % already -1, +1
        errList(ALGO,k,1) = f;
        errList(ALGO,k,2) = pc;
        errList(ALGO,k,3) = pc_test;
        if ~mod( k, 25 )
            fprintf('Iter %3d, train accuracy %.2f%%, test accuracy %.2f%%, objective %.2f\n', ...
                k, pc*100, pc_test*100, f );
        end
    end
end
%%
figure(1); clf;
% subplot(1,3,1)
offset = min( min( errList(:,:,1) ) )-1e-1;
semilogy( errList(1,:,1) - offset, 'linewidth',2 )
hold all
semilogy( errList(2,:,1) - offset, 'linewidth',2 )
semilogy( errList(3,:,1) - offset, 'linewidth',2 )
semilogy( errList(4,:,1) - offset, 'linewidth',2 )
semilogy( errList(5,:,1) - offset, 'linewidth',2 )
semilogy( errList(6,:,1) - offset, 'linewidth',2 )
title('SVM Objective fuction - true value');
xlabel('Iteration');
legend( ALGONAMES )
%% Replot, with corrected x-axis
figure(1); clf;
% subplot(1,3,1)
plotFcn = @loglog;
plotFcn( errList(1,:,1) - offset, 'linewidth',2 )
hold all
plotFcn( linspace(0,maxIts/N,maxIts), errList(2,:,1) - offset, 'linewidth',2 )
plotFcn( linspace(0,maxIts/N,maxIts), errList(3,:,1) - offset, 'linewidth',2 )
plotFcn( linspace(0,maxIts/minibatch_n,maxIts), errList(4,:,1) - offset, 'linewidth',2 )
plotFcn( 1+linspace(0,maxIts/minibatch_n,maxIts), errList(5,:,1) - offset, 'linewidth',2 )
plotFcn( 1:2:(2*maxIts), errList(6,:,1) - offset, 'linewidth',2 )
title('SVM Objective fuction');
xlabel('Epoch');
legend( ALGONAMES )
% xlim([0,3]);

%% Look at miss-classification rate
figure(1); clf;
errMetric   = 2; % train
% errMetric   = 3; % test
% plotFcn = @semilogy;
plotFcn = @plot;
plotFcn( 1-errList(1,:,errMetric), 'linewidth',2 )
hold all
plotFcn( 1-errList(2,:,errMetric), 'linewidth',2 )
plotFcn( 1-errList(3,:,errMetric), 'linewidth',2 )
plotFcn( 1-errList(4,:,errMetric), 'linewidth',2 )
plotFcn( 1-errList(5,:,errMetric), 'linewidth',2 )
plotFcn( 1-errList(6,:,errMetric), 'linewidth',2 )
title('Error, training data');
% title('Error, testing data');
ylabel('Missclassification rate');
legend( ALGONAMES )
xlabel('Iteration');
ylim([0,.15]);
%% Look at miss-classification rate, corrected axis
figure(1); clf;
errMetric   = 3;
plotFcn = @loglog;
plotFcn( 1-errList(1,:,errMetric), 'linewidth',2 )
hold all
plotFcn( linspace(0,maxIts/N,maxIts), 1-errList(2,:,errMetric), 'linewidth',2 )
plotFcn( linspace(0,maxIts/N,maxIts), 1-errList(3,:,errMetric), 'linewidth',2 )
plotFcn( linspace(0,maxIts/minibatch_n,maxIts), 1-errList(4,:,errMetric), 'linewidth',2 )
plotFcn( 1+linspace(0,maxIts/minibatch_n,maxIts), 1-errList(5,:,errMetric), 'linewidth',2 )
plotFcn( 1:2:(2*maxIts), 1-errList(6,:,errMetric), 'linewidth',2 )
title('Error, testing data');
ylabel('Missclassification rate');
legend( ALGONAMES )
xlabel('Epoch');
%% Visualize separating hyperplane
clf;
imagesc( reshape(w(1:end-1),28,28) )
