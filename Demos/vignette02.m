function vignette02

%{
Demo to show the effect of randomized perturbations on the speed
of sorting algorithms

Stephen Becker
%}


rng(0); % make it reproducible

n   = 100;
% x   = randn(n,1);
% x   = (1:n)'; % good for bubble_sort, bad for quick_sort
x   = (n:-1:1)'; % bad for bubble_sort, bad for quick_sort


% == Gather systematic data ==

nList   = round( logspace( 1, 3, 10) );
N       = length( nList );
nReps   = 10;

[bubble, quick]                 = deal( zeros(N,1) );
[bubbleRandom, quickRandom]     = deal( zeros(N,nReps) );
for ni  = 1:N
    n   = nList(ni);
    fprintf('%d of %d trials\n', ni, N );
    
    x   = (n:-1:1)';
    
    a_less_than_b(); % zero out the conter
    y   = bubble_sort(x);
    bubble(ni)  = a_less_than_b(); % number of comparisons
    
    % issorted( y )  % <-- if you don't trust the implementation, check it!
    
    a_less_than_b(); % zero out the conter
    y   = quick_sort(x);
    quick(ni)  = a_less_than_b(); % number of comparisons
    
    for r = 1:nReps
        x   = x(randperm(n));
        
        a_less_than_b(); % zero out the conter
        y   = bubble_sort(x);
        bubbleRandom(ni,r)  = a_less_than_b(); % number of comparisons
        
        a_less_than_b(); % zero out the conter
        y   = quick_sort(x);
        quickRandom(ni,r)  = a_less_than_b(); % number of comparisons
        
    end
    
end
% ==  Plot ==
figure(1); clf;
loglog( nList, bubble, 'o-','linewidth',2,'markersize',8);
hold all
loglog( nList, quick, 's:','linewidth',2,'markersize',8);
loglog( nList, mean(bubbleRandom,2), '*:','linewidth',2,'markersize',8);
loglog( nList, mean(quickRandom,2), 'd:','linewidth',2,'markersize',8);
loglog( nList, nList.^2/2, '--' )
loglog( nList, nList.*log2(nList), '.-' )
set(gca,'fontsize',18);
lh=legend('Bubble Sort','Quicksort','Bubble (randomized)',...
    'Quick (randomized','n^2/2','n log(n)','location','northwest');
xlabel('Length of list');
ylabel('Number of comparisons');
grid on

end % end of main function

%%%%%
% SUBROUTINES
%%%%%

function x = bubble_sort( x )
% x = bubble_sort(x)
%   sorts x in increasing order

n   = length(x);

for iteration = 1:n-1
    
    endIndex    = n - iteration;
    
    alreadySorted   = true;
    
    for i = 1:endIndex
        
        if a_less_than_b( x(i+1), x(i) )
            % swap them:
            tmp = x(i+1);
            x(i+1) = x(i);
            x(i)   = tmp;
            alreadySorted   = false;
        end
        
    end
    if alreadySorted
        break;  % early return
    end
end
end


function x = quick_sort( x )
% x = quick_sort(x)
%   sorts x in increasing order

n   = length(x);

% check for an early return (e.g., the base case in recursion)
%   (a real implementation would have a larger base case)
if n <= 1
    return;
end

% Pick a pivot:
pivot   = x(n); % the last element
% Note: quick sort is never implemented this way, since
%   it uses tricks so that it doesn't have to use extra memory.
smallerList     = [];
largerList      = [];
for i = 1:(n-1)
    value   = x(i);
    if a_less_than_b( value, pivot )
        smallerList(end+1)   = value;
    else
        largerList(end+1)    = value;
    end
end

% Now, recurse
smallerList     = quick_sort( smallerList );
largerList      = quick_sort( largerList );

% and combine
x   = [smallerList; pivot; largerList];
end


function y = a_less_than_b( a, b )
% y = a_less_than_b( a, b )
%   returns "true" if a < b
%   and "false: if a >= b

persistent counter
if nargin == 0
    y   = counter;
    counter = [];
    return;
end
if isempty( counter ), counter = 0; end

% main code:

y   = (a < b );
counter     = counter + 1;
end