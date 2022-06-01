% DiscriminationSim.m
%
% Created 6/00 FMR, Tutorialized JLG 6/00.
%
% In a two alternative forced choice task where humans are
% required to judge whether a very dim flash of light occurs in one
% of two intervals, humans display an amazing ability to
% detect very dim flashes of light. Humans can reliably detect a
% dim light that produces 5-10 photon absorptions spread over 500
% rods.  This ability approaches limits set by noise in the rods.
%
% The rod noise consists of continuous fluctuations and occasional 
% discrete photon-like events.  For just detectable flashes only a 
% small fraction of the rods receive photons, while all generate noise.
% Bipolar cells pool responses from many rods.  How should the bipolar 
% these responses?  Since only a small fraction of the rods inputs to 
% the bipolar receive a photon, but all of the rods are generating noise, 
% identifying those rods that have generated a photon-like response 
% (either due to photon absorption or thermal activation of rhodopsin) 
% could improve the sensitivity of the bipolar cell.  Thus averaging the 
% responses of a number of noisy rod cells with a small number of rods 
% which have absorbed photons could cause the signal to be swamped by 
% continuous noise. How can primates achieve such reliable detection at 
% very low light levels?  This tutorial will show how a static nonlinearity 
% (such as a power-law nonlinearity) applied to the output of each rod
% signal before pooling by a bipolar cell can increase detectability.

% simulation parameters

NumIterations = 100;			% number of simulated trials
SegLength = 100;				% number of time points in each rod response
RodPoolDimension = 12;			% simulated rod array is square with this many rods along one side
MeanPhosA = 0.02;				% mean number of photons absorbed in first response group ('correct' responses)
MeanPhosB = 0.00001;			% mean number of photons absorbed in second response group ('incorrect' responses)
FlashTimeA = 1;					% time point in which flash occurs in first response group
FlashTimeB = 1;					% time point in which flash occurs in second response group
shuffles = 5;					% number of shuffles for bootstrap resampling
Verbose = 1;					% prints lots of crap to command window
RFRadius = 4;					% radius of circularly symmetric Gaussian receptive field
Power = 1;						% parameter of power-law nonlinearity applied to each rod response before pooling

% First we will generate rod responses for our two alternative forced
% choice task. In this case we will create 100 trials where
% we present the stimulus to a 12 x 12 square matrix
% of rod cells. These responses will be stored in a matrix
% which contains all 144 cells response for 100 time points and 100
% iterations. A similar matrix is created for when there
% is no stimulus and contains only the noise inherent in rods.
%
% The simulated rod responses incorporate several
% different components of noise. The first is that if the stimulus has
% say five photons of light, the rod array will on average absorb
% five photons, but the actual number of photons absorbed can be
% modeled as a poisson distributed random variable with mean five.
% A second source of noise is thermal noise. Because of thermal noise
% rhodopsin can randomly change conformation and act as if it had
% just absorbed a photon. This type of noise can also be modeled
% as a poisson distributed variable and gets added to the signal
% created by the actual photon absorptions. There is also noise
% inherent in the signal transdunction pathway that transduces the
% change in rhodopsin confromation to a decrease in cGMP concentrations
% and a decreasing channel open probability. These types of noise
% can be modelled as gaussian distributed with a short temporal
% envelope.
%
% Note that running this code will take quite a bit of time. 
RodResponseSignal = SimulateRodArrayResponse(RodPoolDimension, NumIterations, MeanPhosA, SegLength, FlashTimeA);
RodResponseNoise = SimulateRodArrayResponse(RodPoolDimension, NumIterations, MeanPhosB, SegLength, FlashTimeB);

% Now that we have simulated our response, let's take a look at a single
% neuron's response when there is a stimulus and when there is not a
% stimulus for all 100 of our trials.
%
% What you should see is that this is a complete and utter mess. Over
% 100 presentations of the stimulus, by looking at a single rod there is
% no way to tell the difference between when there was a signal and
% when there was nothing. However, this should not suprise you since
% the probability of this particular rod absorbing a photon is quite
% small, so there may only be a few trials where there was actually
% a response. So, now we have to pool the responses over the whole
% rod array and see what we get.
figure;
for i = 1:NumIterations
  plot(squeeze(RodResponseSignal(6,6,:,i)),'k');hold on
  plot(squeeze(RodResponseNoise(6,6,:,i)),'r');
end
xlabel('time');ylabel('rod response');legend('Rod signal','Rod noise');

% So, now we want to pool the responses of many rods together by
% simply linearly summing the rod responses together.  We want this
% to look somewhat the way a bipolar cell might sum together responses,
% so we assume that our "bipolar" cell will have some sort of receptive
% field which will weigh rod responses over a range of spatial positions.
% We are going to assume our bipolar cell has a gaussian circularly symmetric
% receptive field. To calculate this receptive field we create a synaptic
% weight matrix. 
for XLoc = 1:RodPoolDimension
  for YLoc = 1:RodPoolDimension
	DistanceToRFCenter = (XLoc - RodPoolDimension/2).^2 + (YLoc - RodPoolDimension/2).^2;
	RFWeightMatrix(XLoc, YLoc) = exp(-DistanceToRFCenter / (2 * RFRadius.^2));
  end
end
close;imagesc(RFWeightMatrix);axis off;title('synaptic weights of bipolar cell');

% Now that we have the synaptic weights of our bipolar cell, we
% can pool the responses of all the rods together for our bipolar
% cell and see what our expected bipolar cell response will look like
% on a trial by trial basis. What you should be able to see is that
% summing together the response linearly produces two responses that
% appear quite noisy and very hard to distinguish from each other.
% That is, the noise in the rods is really defining the pooled response
% in the bipolar cell.
figure;subplot(2,1,1)
respSignal(1:NumIterations,1:SegLength) = 0;
respNoise(1:NumIterations,1:SegLength) = 0;
for i = 1:NumIterations
  for XLoc = 1:RodPoolDimension
	for YLoc = 1:RodPoolDimension
  	  temp(1:SegLength) = RodResponseSignal(XLoc, YLoc, 1:SegLength, i);
	  respSignal(i,1:SegLength) = (respSignal(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp)')';
	  temp(1:SegLength) = RodResponseNoise(XLoc, YLoc, 1:SegLength, i);
	  respNoise(i,1:SegLength) = (respNoise(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp)')';
	end
  end
  plot(respSignal(i,:),'k');hold on;plot(respNoise(i,:),'r');
end
title('linear pooling of rod responses, trial by trial basis');
subplot(2,1,2);
stddev = mean(std(respNoise));
plot(mean(respSignal)./stddev,'k','LineWidth',3);hold on;plot(mean(respNoise)./stddev,'r','LineWidth',3);
xlabel('time');ylabel('bipolar response');legend('Signal','Noise');
title('Mean response over 100 presentations scaled to stddev of noise');

% So, our linear pooling of the rod responses seems to be doing
% a poor job of pulling the signal out of the noise. Can we use
% a different pooling mechanism which does a better job? One
% simple idea is to add a static power nonlinearity to the rod
% response. For example, at each time point in each rod response
% we could take a power of the response and then sum together all these
% responses. Will this help make it easier to discriminate
% between the signal and noise?
figure;subplot(2,1,1);
power = 4;
respSignal(1:NumIterations,1:SegLength) = 0;
respNoise(1:NumIterations,1:SegLength) = 0;
for i = 1:NumIterations
  for XLoc = 1:RodPoolDimension
	for YLoc = 1:RodPoolDimension
  	  temp(1:SegLength) = RodResponseSignal(XLoc, YLoc, 1:SegLength, i);
	  respSignal(i,1:SegLength) = (respSignal(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp.^power)')';
	  temp(1:SegLength) = RodResponseNoise(XLoc, YLoc, 1:SegLength, i);
	  respNoise(i,1:SegLength) = (respNoise(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp.^power)')';
	end
  end
  plot(respSignal(i,:),'k');hold on;plot(respNoise(i,:),'r');
end
title('linear pooling of rod responses to the 4th power, trial by trial basis');
subplot(2,1,2);
stddev = mean(std(respNoise));
plot(mean(respSignal)./stddev,'k','LineWidth',3);hold on;plot(mean(respNoise)./stddev,'r','LineWidth',3);
xlabel('time');ylabel('bipolar response');legend('Signal','Noise');
title('Mean response over 100 presentations scaled to stddev of noise');

% What you should see from the previous figure is that adding this
% nonlinearity seems to have brought out the signal in many
% of the signal trials. The reason is that this type of nonlinearity
% stretches larger responses more than smaller responses. Since
% the smaller responses tend to be noise and the larger responses
% tend to be signal, it helps to accentuate the difference between
% signal and noise.

% So exactly how much better did things get? In particular we might
% like to know how this difference in response could be used
% in a psychophysical experiment like the two alternative forced
% choice experiment discussed at the beginning of this tutorial.
% So the experiment is that a stimulus appears in one of two
% intervals and the observer must specify which interval they
% perceived the stimulus. Signal detection theory for one
% dimension says that for the stimulus conditions you get
% a distribution of responses and for the noise conditions
% you get a different distribution of "responses". On a single
% trial, the observer will have one response pulled from
% the signal distribution and another response pulled from
% the noise distribution. By comparing the magnitude of the
% responses and choosing the larger one, the observer can
% decide which interval the signal appeared in. This works
% great when your "response" is one number, say the mean 
% response over time of a bipolar cell. But then we are basically
% throwing out all of the temporal properties of the response.
% How can we generalize the classic one dimensional signal
% detection theory to a stimulus which has many dimensions,
% in this case one dimension for the magnitude of the response
% at each time point? In effect, what we want to do is take
% this multidimensional response and project it on to a one
% dimensional line so that we can again think of it in the
% context of signal detection theory. Of course a simple way
% to do this is just to take the mean response. However this
% weighs each time point the same. Want might be better is
% to take a weighted average, where we weigh points in the
% response where it is likely that there is a large difference
% between stimulus and no stimulus conditions greater then
% when the expected difference is very small. To get good
% weights (under gaussian assumptions), we use the following
% procedure:
%
% We make an n-dimensional space where each axis is the response
% at a different time point. Then we plot the responses for trials
% where the stimulus was present, which we can think of as a gaussian
% blob like cluster of points. And then the responses for trials where
% the stimulus was absent and that should be a separate cluster of
% points. We draw a line between the means of these two blobs. Now
% when we have a new trial, we want to know how far along the
% line that point lies, does it lie closer to the stimulus cluster
% center or the no stimulus cluster center? By taking the dot product
% of the new point with the line, we get the projection of the
% point on to the line, giving us a single number which we can
% now compare with a single number computed in the same way from
% the response in the other interval.
% 
% First we calculate the vector which is the difference between the means.
% Since we want the dot product of this discriminant vector with
% the responses on trials, this is the weight vector which we
% will use to weigh the different time points of the response.
% We calculate this for the last pooling we did where we took
% the 4th power of the rod response.
Discriminant = zeros(1,SegLength);
for i = 1:floor(NumIterations/2)
  Discriminant = Discriminant + respSignal(i,1:SegLength) - respNoise(i,1:SegLength);
end
figure;plot(Discriminant);xlabel('time');ylabel('weight');

% You should see that the discriminant is basically giving high
% weights to the points early in the response and low weights
% to the points late in the trial where the responses are expected
% to be about the same.
%
% Now using this discriminant let's try our 2AFC experiment by
% randomly picking one stimulus present response and one stimulus
% absent response from our calculated distributions. Then we calculate
% what answer we would give using our discriminant. Do this many
% times and then we get an estimate of what the expected percent
% correct would be.
correct = 0;incorrect = 0;
for i = ceil(NumIterations/2+.5):NumIterations
  signalResp = sum(respSignal(ceil(rand*NumIterations),:) .* Discriminant);
  noiseResp = sum(respNoise(ceil(rand*NumIterations),:) .* Discriminant);
  if (signalResp > noiseResp)
	correct = correct + 1;
  else
    incorrect = incorrect + 1;
  end
end
Pcorrect = correct/length(ceil(NumIterations/2+.5):NumIterations)

% So you should see that we get a pretty good probability correct
% for this pooling rule. Now that we have a way of describing
% how good each pooling rule is we go back and calculate the 
% probabiliy correct for different exponents and see how
% good the pooling rules are. This will probably take a bit
% of time to calculate.
power(1:5) = 0;
numbootstraps = 4;
Pcorrect(1:5) = 0;
for n=1:5
  power(n) = 2*n-1;
  % calculate the pooled responses.
  respSignal(1:NumIterations,1:SegLength) = 0;
  respNoise(1:NumIterations,1:SegLength) = 0;
  for i = 1:NumIterations
    for XLoc = 1:RodPoolDimension
	  for YLoc = 1:RodPoolDimension
	    temp(1:SegLength) = RodResponseSignal(XLoc, YLoc, 1:SegLength, i);
	    respSignal(i,1:SegLength) = (respSignal(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp.^power(n))')';
	    temp(1:SegLength) = RodResponseNoise(XLoc, YLoc, 1:SegLength, i);
	    respNoise(i,1:SegLength) = (respNoise(i,1:SegLength)' + RFWeightMatrix(XLoc, YLoc) * (temp.^power(n))')';
	  end
	end
  end

  correct = 0; incorrect = 0;
  for i = 1:numbootstraps
    % calculate discriminant
    shuffled = ShuffleList(1:NumIterations);
	Discriminant = zeros(1,SegLength);
    for i = floor(1:NumIterations/2)
      Discriminant = Discriminant + respSignal(shuffled(i),1:SegLength) - respNoise(shuffled(i),1:SegLength);
    end

    % calculate proability correct
    for i = ceil(NumIterations/2+.5):NumIterations
      signalResp = sum(respSignal(shuffled(i),:) .* Discriminant);
      noiseResp = sum(respNoise(shuffled(i),:) .* Discriminant);
      if (signalResp > noiseResp)
	    correct = correct + 1;
      else
        incorrect = incorrect + 1;
      end
    end
  end
  Pcorrect(n) = correct/(length(ceil((NumIterations/2+.5):NumIterations))*numbootstraps);
end
figure;plot(power, Pcorrect);xlabel('exponent');ylabel('Proability correct');

% What you might notice is that if you increase the exponent too
% high, you also tend to increase the noisiness of the system and
% reduce your detectability, most likely because you create noise
% in your estimation of the discriminant.

% So, what we have found is that have a static exponent nonlineaity
% helps you discriminate pooled rod responses. How might this
% be implemented in neurons? One possibility is that the synapse
% between the rod and the bipolar cell causes this nonlinearity
% through some relationship, say between the amount of calcium
% in the synaptic terminal and the probability of release. Is 
% this what actually happens? Well, the simple experiment to
% do, stimulating a rod and recording from a bipolar cell is
% technically very difficult to do given the small dimensions
% of these cells in the primate retina. So, the answer remains
% to be seen.

