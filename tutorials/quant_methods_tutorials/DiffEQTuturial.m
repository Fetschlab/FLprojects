%****************************************************************************
% Differential equations tutorial for NeuBeh/PBIO 545, Winter 2003
%
% Created 12/02 Fred Rieke
% Revisions:
%
%****************************************************************************


% useful default stuff
set(0, 'DefaultAxesFontName','Palatino')
set(0, 'DefaultAxesFontSize', 12)
colormap([0 0 0])

%****************************************************************************
% Example 1: second-messenger cascades and phototransduction
%	The examples below build up to a model for phototransduction based
% on known chemical interactions.  We start with some of the building
% blocks for a more complicated model and work up to a full model for
% phototransduction.  Many of the examples below are much more generally
% applicable - e.g. the simple feedback models.  Phototransduction is a 
% nice example because we can construct a model based on known mechanisms
% with known properties.
%****************************************************************************

%----------------------------------------------------------------------------
% Consider creation of a substance x by another y.  For example, y could be an active
% receptor and x its downstream effector.  The rate of creation of x 
% is proportional to amount of active y and x decays with a rate 
% constant alpha.  The differential equation describing this situation is:
% 	dx/dt = y - alpha * x
% This does not have a unique solution unless we also specify an initial 
% condition.  We'll use x=0 at time 0.

% The approach we will take to solving this equation is to turn the differential
% equation into a difference equation, which we'll solve for finite time steps.  In
% doing this we will assume that the behavior of x is determined entirely by 
% its first derivative for our finite time step.  This amounts to a Taylor series
% approximation to x(t).  The accuracy of the approximation improves as the time
% step becomes smaller.  If this is unclear come back to it after going through the 
% program below.  

figure(1);
x(1) = 0;				% initial condition
alpha = 20;				% rate constant for decay of x (in 1/sec)
TimeStep = 0.001;		% time step for difference equation (in sec)
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
NumPts = 1000;			% total points to simulate

% initialize y; in this case y is a simple step
y(1:PrePts) = 0;
y(PrePts+1:PrePts + StmPts) = 1;
y(PrePts + StmPts + 1:NumPts) = 0;
% plot y
tme = 1:NumPts;
tme = (tme - PrePts) * TimeStep;
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input y');
title('activation of x by y, no feedback');

% Calculate x by making the differential equation above into a difference equation.  
% Rather than solving with time as a continuous variable, we discretize time into 
% finite steps of length TimeStep.  We will approximate x in each time step.  Thus
% the derivative
%	 dx/dt
% becomes
% 	 [x(n) - x(n-1)] / TimeStep
% where x(n) is the value of x in the nth time bin (i.e. at time n*TimeStep).
% Note that in the limit where TimeStep goes to 0 this is the definition of
% a derivative.  Make sure you see this connection as it is at the core of 
% how differential equations are solved numerically.
% So now our original differential equation becomes
%	[x(n) - x(n-1)] / TimeStep = y(n-1) - alpha * x(n-1)
% or, solving for x(n),
%	x(n) = x(n-1) + TimeStep * [y(n-1) - alpha * x(n-1)]
% This gives us an update rule for x - given x and y in time bin n-1, we 
% can compute x in time bin n.  
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt-1) - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output x');

% Some suggested modifications:
% (1) play with the rate constant alpha and make sure you can explain what happens.
% (2) try different shaped inputs (other than step).
% (3) play with TimeStep to see over what range of time bins the numerical solution
%	  is accurate.

%----------------------------------------------------------------------------
% Let's elaborate this example a little.  What if x also has a rate of spontaneous activation,
% in addition to activation by y?  Now our differential equation becomes
%	dx/dt = y + beta - alpha * x
% where beta is the rate of spontaneous activation.
% We'll solve this in the same way.

figure(2);
beta = 1;				% rate of spontaneous activation of x
x(1) = beta / alpha;	% new initial condition

% plot input y
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input');
title('activation of x by y, no feedback');

% calculate x by making differential equation into difference equation
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt-1) + beta - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output');

% (1) Make sure you understand the form of the difference equation in the loop above.
% (2) What changes in this case?  Why?  
% (3) Can you explain the change from the differential equation?
% (4) Why is beta/alpha is reasonable initial condition?  

%----------------------------------------------------------------------------
% The examples above can both be solved analytically, providing a useful check to 
% the numerical solution.  Now let's add a feedback term, which make the analytical
% solutions difficult.  So now active x will feedback to modify the 'effective' 
% activity of y (that is rate of creation of x by y).  Now the differential
% equation becomes
% 	dx/dt = y * (1+gx)^n - alpha * x,
% where g is a constant describing the gain of the feedback and and n is
% an exponent determining the linear or nonlinear behavior of the feedback signal.
% There are certainly other ways a feedback mechanism could work (and correspondingly
% different differential equations), but this is one reasonable form.

figure(3);
x(1) = 0;						% initial condition
Power = -1;						% feedback power
g = 20;							% feedback gain

% plot input y
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input');
title('activation of x by y, feedback');

% solve difference equation
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt) * (1 + g * x(pnt-1))^Power - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output');

% (1) How does the behavior compare to the case without feedback?
%	  Compare both the amplitude and kinetics of x.
% (2) How does effect of feedback change as power changed?  Why? 
% (3) How does changing g change things?  Why?
% (4) What happens when you change TimeStep?  Why?

%----------------------------------------------------------------------------
% Now let's consider combinations of a couple steps.  Anticipating building up
% to a model for phototransduction we'll change from the generic variables x and y
% to r (rhodopsin activity) and p (phosphodiesterase activity).  Light produces
% an electrical signal in the photoreceptor by activating rhodopsin which then
% activates phosphodiesterase (through transducin).  Active phosphodiesterase 
% breaks down cGMP and reduces the membrane current through cGMP-gated channels.  
% Start by considering phosphodiesterase activation.  First, we'll assume that the 
% rhodopsin activity (its ability to activate phosphodiesterase through transducin) 
% obeys
%	dr/dt = -sigma*r
% where sigma is the rate constant for the decay of rhodopsin.  Then 
% we'll assume that the phosphodiesterase activity obeys
%	dp/dt = r + eta - phi * p 
% where eta represents spontaneous phosphodiesterase activation and
% phi is the rate constant for decay.

figure(1)
sigma = 5;				% rhodopsin activity decay rate constant (1/sec)
phi = 5;				% phosphodiesterase activity decay rate constant (1/sec)
eta = 10;				% phosphodiesterase activation rate constant (1/sec)
r(1) = 1;				% initial condition for r
p(1) = eta/phi;			% initial condition for p

NumPts = 1000;			% number of points to simulate
TimeStep = 0.001;		% time step 

tme = 1:NumPts;
tme = tme * TimeStep;

% solve difference equation
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
end

% plot time course of rhodopsin activity
subplot(2, 1, 1);
plot(tme, r);
xlabel('time (sec)')
ylabel('rhodopsin activity')
% plot time course of phosphodiesterase activity
subplot(2, 1, 2);
plot(tme, p);
xlabel('time (sec)')
ylabel('pde activity')

% (1) why do we choose an exponential decay for the shape of rhodopsin's activity?  What might
%	  change that?
% (2) explore different combinations of sigma and phi and their impact.

%----------------------------------------------------------------------------
% The code above describes the 'activation' part of phototransduction - i.e. how
% light activation of rhodopsin leads to activation of phosphodiesterase.  Now 
% lets add in the steps linking that to a change in current, and the steps 
% required for the light response to recover.  The role of activated phosphodiesterase
% is to hydrolyze cGMP.  Another enzyme, guanylate cyclase, synthesizes cGMP.  So the 
% cGMP concentration g is controlled by a balance of synthesis (at a rate s) and hydrolysis
% (at a rate pg):
%	dg/dt = s - pg
% The membrane current depends on the third power of the cGMP concentration
%	I = k g^3
% This is all we need for a simple phototransduction model.  Below we add a feedback
% that controls the rate of synthesis.

figure(1)
gdark = 15;				% concentration of cGMP in darkness
cgmp2cur = 8e-3;		% constant relating cGMP to current

tme = 1:NumPts;
tme = tme * TimeStep;
g(1) = gdark;							% initial condition
s(1:NumPts) = gdark * eta/phi;		% steady state: synthesis = hydrolysis
p(1) = eta/phi;

% solve difference equations for each component
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
end

% determine current change
cur = cgmp2cur * g.^3;

% plot current, pde activity, synthesis rate, cGMP concentration
subplot(4, 1, 1);
plot(tme, cur);
xlabel('time (sec)');
ylabel('current');
subplot(4, 1, 2);
plot(tme, p);
xlabel('time (sec)');
ylabel('pde activity');
subplot(4, 1, 3);
plot(tme, s);
xlabel('time (sec)');
ylabel('synthesis rate');
subplot(4, 1, 4);
plot(tme, g);
xlabel('time (sec)');
ylabel('[cGMP]');

%----------------------------------------------------------------------------
% Now let's add a calcium feedback to the model above.  Calcium enters the outer
% segment through the cGMP-gated channels.  Hence the rate of calcium influx
% is proportional to the current flowing.  Calcium is removed by a Na+/K+,Ca2+ 
% exchanger.  The rate of removal is proportional to the calcium concentration:
%	dc/dt = qI - beta c
% where c is the calcium concentration, q is the proportionality constant between
% changes in calcium and the current, and beta is the rate constant for calcium
% removal.  
% Calcium acts on the rate of cGMP synthesis:
%   s = smax / (1 + (c/K)^h)
% where smax is the maximum rate and K and H are constant (affinity and cooperativity)
% describing the strength of the feedback.  Otherwise the model is the same.

figure(2)
cdark = 0.5;			% dark calcium concentration
beta = 20;				% rate constant for calcium removal in 1/sec
hillcoef = 4;			% cooperativity
hillaffinity = 0.3;		% affinity

cur2ca = beta * cdark / (cgmp2cur * gdark^3);				% get q using steady state
smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);		% get smax using steady state

tme = 1:NumPts;
tme = tme * TimeStep;

% initial conditions
g(1) = gdark;
s(1) = gdark * eta/phi;		
c(1) = cdark;
p(1) = eta/phi;

% solve difference equations
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
	s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
end
% determine current change
cur = cgmp2cur * g.^3;

% plot current, pde, synthesis, cGMP and calcium
subplot(5, 1, 1);
plot(tme, cur);
xlabel('time (sec)');
ylabel('current');
subplot(5, 1, 2);
plot(tme, p);
xlabel('time (sec)');
ylabel('pde activity');
subplot(5, 1, 3);
plot(tme, s);
xlabel('time (sec)');
ylabel('synthesis rate');
subplot(5, 1, 4);
plot(tme, g);
xlabel('time (sec)');
ylabel('[cGMP]');
subplot(5, 1, 5) 
plot(tme, c)
xlabel('time (sec)');
ylabel('[calcium]');

% (1) make sure you understand the relation between the steady state conditions
%	  and the constants q and smax.
% (2) play with the various parameters and see how the alter the calculated 
%	  light response.  Can you explain why things change they way they do?
% (3) the model will generate damped oscillations for some values of beta.  Why?

%****************************************************************************
% Example 2: two-state systems and channel gating particles
%****************************************************************************

% two-state system - e.g. HH gating particle
%	dp/dt = alpha n - beta p
%	p+n = 1

n(1) = 1;
alpha = 50;
beta = 50;
NumPts = 1000;
TimeStep = 0.0001;

p(1:NumPts) = 0;
for pnt = 2:NumPts
	p(pnt) = p(pnt-1) + TimeStep * (n(pnt-1) * alpha - p(pnt-1) * beta);
	n(pnt) = 1 - p(pnt);
end

tme = 1:NumPts;
tme = tme * TimeStep;
subplot(2, 1, 1);
plot(tme, p);
subplot(2, 1, 2);
plot(tme, n);

% (1) explain ss dependence on alpha and beta
% (2) add inactive state entered from permissive state
