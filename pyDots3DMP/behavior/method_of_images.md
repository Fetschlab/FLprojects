# Method of Images Approach and Anti-Correlated Accumulators

### First-passage Monte Carlo simulation
For a 2-D accumulator, we can empirically simulate the decision variable 
using accumulation of 1-D gaussians - (this is standard first-passage Monte Carlo), 
then take the identity and time of first bound crossing.
Repeating this many times and averaging would give a (relatively) unbiased estimate of the probability

To prepare for MOI, imagine representing the two accumulators not on a single axis of 
[-Inf->bound] vs time, but a 2-D cartesian grid, where the diffusion particle starts at 
[-bound, -bound]. The [-Inf->0] segments of x and y-axes now represent the two boundaries.

### Direct solution
Now, we want to go after the probability distributions more directly,
so that we don't have to simulate N trials to approximate them. This means we 
want to enforce absorbing boundaries i.e. density of the diffusion particle is zero 
at ALL points at the boundaries.
Absorbing boundaries cannot be enforced directly with pdf calculations, 
the method of images provides a way around this - images are reflected copies of the
free-space solution that cancel out/negate the probability that would "leak" across the boundary.

How do we determine the number of images needed?

For a 2-D cartesian grid with 4 quadrants and the lower left being our operating space,
there are 3 places probability could leak:
    1) upper left (-x, +y)
    2) upper right (+x, +y)
    3) lower right (+x, -y)
which necessitates 3 images (i.e. 2 reflections - one across the y axis, one across x)

HOWEVER, this is only holds when rho = 0 i.e. the two accumulators are orthogonal,
because only then is the diffusion operator invariant under reflection across the x/y boundaries.
When we introduce (anti-)correlation, we need to change basis to the eigen-directions of rho.
Once in this situation, the images do not "close" after only 2 reflections.
The general equation governing how many reflections (k) are needed is:
    rho = -cos(pi/k)
and the number of images = 2k-1

For an anti-correlation of 1/sqrt(2), we therefore need k=4 reflections and thus 7 images.


### To come
More explanation of cdf, pdf calculations (vectorized and loop versions)
