% ----------------------------------------------------------------------
% whiteNoise tutorial
% ----------------------------------------------------------------------
% This is a script which simulates a white noise experiment.  
% In brief, we're going to simulate the neural response to a
% one dimensional visual stimulus (luminance with respect to the
% mean, or contrast) and then reconstruct the mechanism of the
% neural response using the techniques described in the handout
% by E.J. Chichilnisky.  The neural response is simulated as a
% non-homogenous Poisson process whose rate parameter (which will be called
% "nonlinearResp" below) is a linear function of the stimulus (which
% will be called "linearResp" below) put through a static, single-valued
% non-linearity.
% ----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Initialize.
%
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We're going to run the experiment for 250 seconds or slightly more than 
% 4 minutes.  The method works better with longer times, but this works 
% pretty well and it doesn't crash the machine.
duration=250000;            % Time is measured in msecs
stimSampleTime=10;			% 100 Hz monitor

% Our "monitor" (device which generates the stimulus)
% refreshes the display at 100 HZ, or once every 10 msecs.
% At each frame, the intensity of the uniform screen changes:
% it is drawn randomly from a particular probability distribution.
% The exact nature of the probability distribution is not terribly important.
% However, it is important that 1) the distribution be symmetric about
% zero and that 2) the contrast on each frame is uncorrelated with the others
% (this is the definition of a "white" random process).
% As an aside: For stimuli that vary in space and time, 
% the distribution must be symmetric about the origin.  
% For practical purposes, this means it should be Gaussian.

x=clip(1/3*randn(duration/stimSampleTime,1),-1,1);
stimulus=zeros(duration,1);
for i=1:length(x)
  t=10*(i-1)+1;
  stimulus(t:t+9)=x(i)*ones(10,1);
end

plot(stimulus(1:1000))
xlabel('Time (msec)')
ylabel('Stimulus contrast')
set(gcf,'Units','Normalized','Position',[.02 .7 .3 .3])
title('First 1 second of stimulus');

% Compute the linear response using a 3 stage cascade of exponential
% lowpass filters.  You could substitute other linear filters here if you
% wanted to, but this one is pretty simple.  The 3rd row of y is the linear
% response (refer to diffEqTutorial).  This takes couple of minutes to calculate
linearResp=zeros(length(stimulus),3);
tau=15;
for i=1:length(stimulus)-1
  linearResp(i+1,1) = linearResp(i,1) + (1/tau) * (stimulus(i) - linearResp(i,1));
  linearResp(i+1,2) = linearResp(i,2) + (1/tau) * (linearResp(i,1) - linearResp(i,2));
  linearResp(i+1,3) = linearResp(i,3) + (1/tau) * (linearResp(i,2) - linearResp(i,3));
end
% Getting rid of the first- and second-order filtered signals, we only
% want the third one.
linearResp = linearResp(:,3);

% The linear response is just a lowpass filtered version of the stimulus.  The top panel
% of the figure shows the first second of the stimulus, the small middle panel shows
% the impulse response of the linear filter, and the third panel shows the first second
% of the linear response.
set(gcf,'Units','Normalized','Position',[.02 .7 .3 .6])
subplot(2,1,1);
plot(stimulus(1:1000))
title('Stimulus');
set(gca,'XTIck',[])
subplot(2,1,2);
plot(linearResp(1:1000))
title('Linear response');
xlabel('Time (ms)');
axes('position',[.13 .47 .2 .1])
times=[0:200]';
impulseResponse = times.^2 .* exp(-times/tau);
plot(times,impulseResponse/max(impulseResponse),'m-')
set(gca,'XTick',[],'YTick',[],'Box','on','Units','Normalized')
h = text(40,.6,'Imp. Resp.');
set(h,'FontSize',9)

% Under the simple non-linear (SNL) model, the firing rate of a neuron is
% a single-valued non-linear function of an underlying linear response.
% We can pick any such function we want, and, for this example, we've decided
% on the cumulative Gaussian function (the integral of a Gaussian probability
% density function).  
% The way we've parameterized it, this function has three arguments:
% the slope (alpha), the point of inflection (-beta/alpha), and the upper
% asymptote (gamma). This function relates the probability of firing to the 
% linear response.
alpha=25;
beta=-2;
gamma=.15;
subplot(1,1,1)
plot([-.25:.01:.3],gamma*normcdf(alpha*[-.25:.01:.3]+beta,0,1))
set(gca,'XLim',[-.25 .3]); set(gcf,'Units','Normalized','Position',[.02 .7 .3 .3])
ylabel('Probability of firing'); xlabel('Linear response');

% Here we're applying this non-linear transformation on the linear response of
% our simulated neuron.  This takes a few seconds and if you're going to run
% out of memory it's going to happen now.
nonlinearResp=gamma*normcdf(alpha*linearResp+beta);

% We can use this non-linear response to simulate a Poisson-ish spike train...
xr=rand(size(nonlinearResp));
neuralResponse = nonlinearResp>xr;

% ...and we can count up the number of spikes fired by the neuron in each 10 msec
% wide bin (each screen refresh).
spikeCounts=zeros(length(stimulus)/10,1);
for i=1:length(stimulus)/10;
  t=10*(i-1)+1;
  spikeCounts(i)=sum(neuralResponse(t:t+9));
end

% So far, we constructed a white noise stimulus, we linearly filtered it, put this
% linearly filtered signal through a non-linear function to calculate an
% underlying firing rate of the cell, and used this underlying firing rate to 
% simulate spikes coming out of the cell.  Here's the first second of each of these
% functions 
h = figure
set(gcf,'Units','Normalized','Position',[.12 .09 .7 .8])
h = subplot(4,1,1)
plot(linearResp(1:1000),'b-')
title('Linear response')
set(h,'XTick',[],'XTickLabel',[])
h = subplot(4,1,2)
plot(nonlinearResp(1:1000),'r-')
title('Firing probability (nonlinear function of linear response)')
set(h,'XTick',[],'XTickLabel',[])
h = subplot(4,1,3)
plot(neuralResponse(1:1000))
set(gca,'Ylim',[0 2]);
title('# of Spikes (1 ms bins)');
set(h,'XTick',[],'XTickLabel',[])
h = subplot(4,1,4)
plot([1:10:1000],spikeCounts(1:100),'m-')
title('# of Spikes (10 ms bins)');
xlabel('Time (ms)');

% Now we compute the spike-triggered average stimulus.  This is accomplished
% by taking the 30 milliseconds of stimulus immediately preceding each spike
% and adding them together.  This sum is then divided by the total number 
% of spikes fired over the course of the entire experiment to determine the 
% average stimulus preceding a spike.
% This spike-triggered average is, in a sense, a template for what the neuron
% is "looking for" in the stimulus.

totalCount=sum(spikeCounts)
windowSize=30;
a=zeros(windowSize,1);
for i=windowSize+1:length(spikeCounts)
  a = a + spikeCounts(i)*stimulus((i-windowSize)*10:10:i*10-1);
end
a=a/totalCount;
figure
plot(a)

% A powerful result from the theory of white noise analysis states that 
% if  the neuron satisfies the assumptions of the SNL model, and the 
% stimulus is drawn according to a distribution which is symmetric about the origin,
% the spike-triggered average converges, as the time of the experiment goes to 
% infinity, to to the (time-reversed) impulse response of the linear part of the SNL
% system (up to an arbitrary scale factor).  To the extent that neurons are well-
% modeled as SNL systems, this means that we can easily measure the linear
% component.  Because this is a tutorial, we *know* exactly what filtering 
% was done on the stimulus to get the linear response ("linearResp").  Recall
% that it was a cascade of three first-order exponential filters.  Below, we 
% compare the spike-triggered average to the impulse response of the filter we used.
times=[0:400]';
impulseResponse = times.^2 .* exp(-times/tau);
impulseResponse = impulseResponse/max(impulseResponse);
aFlipped=flipud(a);
sImpResponse = impulseResponse(10:10:10*windowSize+9);
clf;
plot([aFlipped,sImpResponse])

% There is a discrepancy between the impulse response of the filter and the 
% (time-flipped) spike-triggered average.  One is a scaled version of the other.
% There is an ambiguity between the magnitude of the linear response and the
% scale of the non-linear function (intuitively, we could get identical neural
% responses either by scaling the linear response or by scaling the subsequent
% non-linear function).

% Here they are again, this time scaled to have the same energy
plot([aFlipped/sum(aFlipped),sImpResponse/sum(sImpResponse)])

% Now we can estimate the non-linear function which relates the linear response
% to the probability of firing a spike.  This is accomplished by plotting the 
% estimated linear response versus the spike count in each 10 ms bin.  It turns 
% out that a simple scatter diagram is pretty uninformative because most of the 
% points overlap each other...

linearEst=zeros(size(spikeCounts));
for i=windowSize+1:length(spikeCounts)
  linearEst(i) = a' * stimulus((i-windowSize)*10:10:i*10-1);
end
plot(linearEst,spikeCounts,'o')

% We can clean this plot up by plotting the *average* number of spikes fired 
% in response to similar linear responses.  First we decide on linear
% response bins...
x=[-.2:.05:.3]; 
% ...and then calculate the averages.
N=zeros(size(x));
se=zeros(size(x));
for i=1:length(x)
  L = logical((linearEst>x(i)) & (linearEst<(x(i)+.05)));
  N(i)=mean(spikeCounts(L));
  se(i)=sqrt(var(spikeCounts((linearEst>x(i)) & (linearEst<(x(i)+.05))))/sum(L));
end
x=x+.025;
errorbar(x,N,se,'ro')
set(gca,'xLim',[min(x), max(x)]);
xlabel('Linear response component');
ylabel('Mean spike count');

% Now we can superimpose the non-linear function that we actually used to 
% determine the spike firing probabilities.  The plot of linear response versus
% mean spike count should have the same shape as this function, but remember, there
% was an arbitrary scale factor relating these two quantities.  Below, we estimate
% this scale factor using least-squares.
xx = normcdf(alpha*x+beta);
gamma = regress(N',xx')
x=[-.4:.005:.4]; 
Nth=gamma*normcdf(alpha*x+beta);
hold on
plot(x,Nth)
hold off