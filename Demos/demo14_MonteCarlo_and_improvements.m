%{
Discusses Monte Carlo in the context of integration:

- There are many ways to integrate functions

- Deterministic "quadrature" rules are fancy Riemann Sums, and
    will work *very well* if the integrand is smooth
    They break down when the integrand is highly oscillatory,
    and/or for high-dimensional integrals. Special versions targeted
    for oscillatory integrals is the subject of current applied math
    research.

- Monte Carlo integration interprets the integral as an expectation
    of a random variable, and draws samples to approximate the true mean
    with a sample mean.
    For a smooth function, Monte Carlo integration is a bad idea because
    classical quadrature rules are much, much better

- Monte Carlo is slow/inaccurate, but the inaccuracy is independent
    of the dimension of the integral. So for large enough dimensions, 
    it makes sense (while in large dimensions, making a deterministic
    grid is impossible since it will be too large)

- Since Monte Carlo is useful sometimes, there are many known techniques
    to make it better. We examine two:
        -- Quasi Monte Carlo, which uses low-discrepancy sequences, and
            inherits some of the advantages and disadvantages from 
            both Monte Carlo and grid/quadrature methods.
            Refs:
            - https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Construction_of_low-discrepancy_sequences
            - "High-dimensional integration: The quasi-Monte Carlo way" by Dick, Kuo
            and Sloan (Acta Numerica, 2013)
        -- Control variates as a means of variance reduction
            Refs: 
            - https://en.wikipedia.org/wiki/Control_variates

Stephen Becker, University of Colorado, April 2019
%}

%% Integrate sin(x)/x from 0 to 1 (e.g. Si(1), Si is Sine Integral)
%{
    The sine integral, Si(z), is the integral of sin(x)/x from 0 to z
    where we define sin(0)/0 to be 0 (consistent with the limit)

    This integral is not known in closed form
    See https://en.wikipedia.org/wiki/Trigonometric_integral#Sine_integral

    How can we approximate it? There are specialized techniques that are
    faster and more accurate than what we will discuss here, but we'll
    treat it via the integral definition and try to numerically
    evaluate the integral.
%}
si  = sinint(1);       % get fairly accurate answer using Matlab's symbolic toolbox
f   = @(x) sinc(x/pi); % equivalent to sin(x)/x and f(0)=0
N   = 1e2+1; % keep it odd for my composite Simpson's to work
xgrid   = linspace(0,1,N);
dx      = xgrid(2)-xgrid(1);
fx      = f(xgrid);
composite_mid   = dx*sum(f(xgrid(2:end)-dx/2)); % open formula
composite_trap  = dx*( sum(fx) -fx(1)/2 - fx(end)/2 );
composite_simp  = dx/3*( fx(1)+fx(end)+ 4*sum(fx(2:2:end-1)) + 2*sum(fx(3:2:end-1)) );
si - composite_mid
si - composite_trap
si - composite_simp

%% 2a visualize discrepancy of random numbers on [0,1]

N   = 1e3;
setA    = sort(rand(N,1));
setB    = [.5*setA(1:2:end); .5 + .5*setA(2:2:end)];

figure(1); clf;
plot( setA, 'linewidth',2 ); hold all; plot( setB, 'linewidth',2 ); 
legend('uniform random','lower discrepancy');
line([0,N],[0,1],'linestyle','--','color','k');
%% more plots
figure(1); clf;
area( smooth(setA - linspace(0,1,N)') ); hold all
ar=area( smooth(setB - linspace(0,1,N)') ); line([0,N],[0,0],'color','k');
ar.FaceAlpha = 0.5; ar.FaceColor = 'r';
legend('uniform random','lower discrepancy');
%% ore plots
clf;
histogram( diff(setA) ); hold all
histogram( diff(setB) ); legend('uniform random','lower discrepancy');
title('Separation distances in random "grid"')
%% Try Monte Carlo evaluation of Si(1)
N   = 1e2;
setA    = sort(rand(N,1));
setB    = [.5*setA(1:2:end); .5 + .5*setA(2:2:end)];

int_MonteCarlo  = mean(f(setA));
int_QuasiMonteCarlo     = mean( f(setB) );

% Add in control variate
% Use sin(x)/x ~ 1 - x^2/6 (first part of Taylor series)
g   = @(x) 1 - x.^2/6;
% The integral (or mean/expectation) of g over [0,1] is:
int_g   = 17/18;
% si - int_g  % already a fairly good approximation
fx      = f(setA);
gx      = g(setA);
% Estimate covariance and variance of gx
cv      = cov(fx,gx);
c       = -cv(1,2)/cv(2,2); % estimated optimal control variate parameter
int_ControlVariate  = int_MonteCarlo + c*(mean(gx)-int_g);
int_ControlVariate_quasi  = int_QuasiMonteCarlo + c*(mean(g(setB))-int_g);
fprintf('\nError is %10.3e for plain Monte Carlo\n',  si - int_MonteCarlo );
fprintf('Error is %10.3e for Quasi Monte Carlo\n', si - int_QuasiMonteCarlo );
fprintf('Error is %10.3e for 2nd order Taylor Series\n', si - int_g );
fprintf('Error is %10.3e for Control-Variate Monte Carlo\n', si - int_ControlVariate );
fprintf('Error is %10.3e for Control-Variate Quasi Monte Carlo\n', si - int_ControlVariate_quasi );
fprintf('Error is %10.3e for quadrature (composite Trapezoidal Rule)\n',si - composite_trap);
fprintf('Error is %10.3e for quadrature (composite Simposon''s Rule)\n',si - composite_simp);