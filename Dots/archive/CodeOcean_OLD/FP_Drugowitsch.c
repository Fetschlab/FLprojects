/*
* arrayProduct.c - example in MATLAB External Interfaces
*
* Multiplies an input scalar (multiplier)
* times a 1xN matrix (inMatrix)
* and outputs a 1xN matrix (outMatrix)
*
* The calling syntax is:
*
*		outMatrix = arrayProduct(multiplier, inMatrix)
*
* This is a MEX file for MATLAB.
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

//Matlab requirment
#include "mex.h"

double pi = 3.14159265358979323846;
const int N = 2;
//int xmeshSize = 301 * 301; // temporary, make this an input later

// You have to write out functions for all simple stuff, like finding determinant, multplying matrix, and so forth
double determinantOfMatrix(double** matrix){
	return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
}

// Function for scalar multiplication expansion 
double** dotMultiplication(double** matrix, double value){
	double** newMatrix = (double **)mxCalloc(2, sizeof(double));
	for (int i = 0; i < N; ++i){
		newMatrix[i] = (double *)mxCalloc(2, sizeof(double));
		for (int z = 0; z < N; ++z){
			newMatrix[i][z] = matrix[i][z] * value;
		}
	}
	return newMatrix;
}
double* dotMultiplication1D(double* vector, double value, int xmeshSize){
	for (int i = 0; i < xmeshSize; i++){
		vector[i] = vector[i] * value;
	}
	return vector;
}
//addVectors
double* addVectors(double* vector1, double* vector2, int xmeshSize){
	//double* newVector = (double *)mxCalloc(xmeshSize, sizeof(double));
	//double* newVector = new double[xmeshSize];
	for (int i = 0; i < xmeshSize; i++){
		vector1[i] = vector1[i] + vector2[i];
	}
	return vector1;
}
// Solving for inside the exponential in the gaussian
double insideExp(double X[2], double* S, double* mu, double** sigma, double t){
	// First find the inv or Sigma
	double det = determinantOfMatrix(sigma);
	double tl = sigma[1][1]/det; //[0][0]
	double tr = -sigma[0][1]/det; // [0][1]
	double bl = -sigma[1][0]/det; // [1][0]
	double br = sigma[0][0]/det; // [1][1]

	// Now do the dirty work of matrix multiplication
	double tempOne = (X[0] - S[0] - mu[0] * t) *tl + (X[1] - S[1] - mu[1] * t) *bl;
	double tempTwo = (X[0] - S[0] - mu[0] * t) * tr + (X[1] - S[1] - mu[1] * t) *br;

	double insideExponential = tempOne * (X[0] - S[0] - mu[0] * t) + tempTwo * (X[1] - S[1] - mu[1] * t);
	return insideExponential;
}

// Make 2D Gaussian Profile
double* twoDGauss(double* X, double* S, double* mu, double** sigma, double t, int xmeshSize){
	double leftHalf;
	// Calculate the Left Side of the Gaussian
	//double** newSigma = dotMultiplication(sigma, t); // Calculate new sigmaT, not needed
	leftHalf = 1 / (2 * pi * sqrt(t*t * determinantOfMatrix(sigma))); // add the (smallest number)^2 possible to this equation. I think you multiply the t with the Sigma, but it comes out since you are doing the determinant

	//Calculate the right side of the gaussian
	double* probabilityVector = (double *)mxCalloc(xmeshSize, sizeof(double));
	// loop through the xmesh grid length (its not 1000 here, let it be flexible; mayhbe it let it be an input)
	for (int j = 0; j < xmeshSize; ++j){
		double tempX[2] = { X[j], X[j + xmeshSize] }; // oh snap this should work!!!
		double tempExp = insideExp(tempX, S, mu, sigma, t); // send the pointer X
		probabilityVector[j] = leftHalf * exp(-1/(2*t) * tempExp); // You might need to fix the T here so it correctly multiplies with the sigma
	}
	//if (newSigma){ mxFree(newSigma); } // Free Data for performance
	return probabilityVector;
}
// Sj_even
double* js_even(double iteration, double alpha, double k, double* s0){
	double* newS = (double *)mxCalloc(2, sizeof(double));
	newS[0] = (1 / sin(pi / k) * sin(iteration*alpha + pi / k) * s0[0]) + (1 / sin(pi / k) * sin(iteration*alpha) * s0[1]);
	newS[1] = (1 / sin(pi / k) * -sin(iteration*alpha) * s0[0]) + (1 / sin(pi / k) * -sin(iteration*alpha - pi / k) * s0[1]);

	return newS;
}
//Sj_Odd
double* js_odd(double iteration, double alpha, double k, double* s0){
	double* newS = (double *)mxCalloc(2, sizeof(double));
	newS[0] = (1 / sin(pi / k) * sin(iteration*alpha) * s0[0]) + (1 / sin(pi / k) * sin(iteration*alpha - pi / k) * s0[1]);
	newS[1] = (1 / sin(pi / k) * -sin(iteration*alpha + pi / k) * s0[0]) + (1 / sin(pi / k) * -sin(iteration*alpha) * s0[1]);

	return newS;
}
// Inside exponential for Weight function
double insideWeightExp(double* mu, double** sigma, double* sj, double* s){
	// First find the inv or Sigma
	double det = determinantOfMatrix(sigma);
	double tl = sigma[1][1] / det; //[0][0]
	double tr = -sigma[0][1] / det; // [0][1]
	double bl = -sigma[1][0] / det; // [1][0]
	double br = sigma[0][0] / det; // [1][1]

	// Now do the dirty work of matrix multiplication (you can probability remove the dynamic part here!!!!)
	//double subS[2] = { sj[0] - s[0], sj[1] - s[1] };
	//double* subS = (double *)mxCalloc(2, sizeof(double));
	//subS[0] = sj[0] - s[0]; subS[1] = sj[1] - s[1];
	// Now make new variable
	//double insideWeightTerm = (mu[0] *tl + mu[1] *bl) * subS[0] + (mu[0] *tr + mu[1] *br) * subS[1];
	double insideWeightTerm = (mu[0] * tl + mu[1] * bl) * (sj[0] - s[0]) + (mu[0] * tr + mu[1] * br) * (sj[1] - s[1]);
	//if (subS){ mxFree(subS); }
	return insideWeightTerm;
}

//Weight values
double weightj(double iteration, double* mu, double** sigma, double* sj, double* s){
	double weight;
	double tempExp = insideWeightExp(mu, sigma, sj, s);
	weight = pow(-1, iteration) * exp(tempExp);
	return weight;
}


/* Main Algorithm -> Right from Matlab, just copied and paste (there are other algorithms I use, keep that in mind) */
long drugo(double* subInd, double delta, int xmeshSize, double* mu, double k, double totalTime, double dt, double* s0, double alpha, double** globalSigma, double urgencyMax, double urgencyTauHalf, double* xytFPFinal){
	/* Variables needed for algorithm*/
	// Variable for Covariance matrix
	//double delta = .01;

	// Parameters that can be dump at the end 
	double *restAdd;
	double *addRest;
	double* evenS;
	double* oddS;
	double weightForT;
	int z;
	double kInd;
	double t;
	int r;
	double* additionalMu = (double *)mxCalloc(2, sizeof(double));

	// Main Drugowitsch Algorithm (You might have to construct other functions within algoritm to break it down easier)
	for (t = 0, r = 0; t < totalTime; t = t + dt, ++r){
		// Here we will calculate the addition of the urgency signal 
		// I have also added 'additionalMu' to every 'mu'
		if (urgencyMax != 0 && urgencyTauHalf != 0){
			additionalMu[0] = mu[0] + urgencyMax * (1 / (t + urgencyTauHalf)); // + (urgencyMax * (t / (t + urgencyTauHalf)))/t
			additionalMu[1]  = mu[1] + urgencyMax * (1 / (t + urgencyTauHalf));
		}
		else if(urgencyMax != 0 && urgencyTauHalf == 0){
            // No need to multiply the Urgency calculation by 'r' because this it is in fact multiply by an 'r' in the 'twoGauss' calculation
            // aY = a(mu + Um) -> aY = a*mu + a*Um... 'a*mu' is the original calculation and linear urgency signal would be 'a*Um'
            additionalMu[0] = mu[0] + (urgencyMax); // * .001); // should this be -
            additionalMu[1] = mu[1] + (urgencyMax); // * .001);
		}
		else{
			additionalMu[0] = mu[0] + 0;
			additionalMu[1] = mu[1] + 0;
		}
        
        addRest = twoDGauss(subInd, s0, additionalMu, globalSigma, t, xmeshSize); // Erase 'addRest' at the end of the loop
        
        // You were mistakenly writing k*2 - 1, which is correctly writing but in C code the for loop requires this to be k*2 ignoring the -1.
		for (kInd = 1; kInd < k * 2; ++kInd){
			// is it even???
			if (round(kInd / 2) == kInd / 2) {
				evenS = js_even(kInd, alpha, k, s0);
				restAdd = twoDGauss(subInd, evenS, additionalMu, globalSigma, t, xmeshSize);
				weightForT = weightj(kInd, additionalMu, globalSigma, evenS, s0);
				if (evenS){ mxFree(evenS); } // Free Data for performance
			}
			// if not than it must be odd
			else {
				oddS = js_odd(kInd, alpha, k, s0);
				restAdd = twoDGauss(subInd, oddS, additionalMu, globalSigma, t, xmeshSize);
				weightForT = weightj(kInd, additionalMu, globalSigma, oddS, s0);
				if (oddS){ mxFree(oddS); } // Free Data for performance
			}
			// Missing one thing here
			// 1) Add the addRest = addRest + weightForT .* restAdd; (done)
			restAdd = dotMultiplication1D(restAdd, weightForT, xmeshSize);
			addRest = addVectors(addRest, restAdd, xmeshSize);
			if (restAdd){ mxFree(restAdd); } // Free Data for performance
		}

		// You need two things 1) Multiply by meshGrid^2 and 2) Reshape... however for reshape you might be better off doing this in Matlab!!!!
		// WAIT YOU HAVE TO ALSO CONNECT EACH TIME, SO BASICALLY STICH THINGS TOGETHER TO FORM A 3D MAT, FIND A SMART WAY TO MAKE IT 2D AND CHANGE IT IN MATLAB
		addRest = dotMultiplication1D(addRest, delta*delta, xmeshSize);

		//mexPrintf("you got here \n");
		//mexEvalString("drawnow;");
		/* MAMA MIA; this is where saving a 3D matrix into a 1D pointer gets a bit difficult*/
		for (z = 0; z < xmeshSize; ++z)
		{
			xytFPFinal[r * xmeshSize + z] = addRest[z];
		}
		//mexPrintf("you're out of xytFinal \n");
		//mexEvalString("drawnow;");
		if (addRest){ mxFree(addRest); } // Free Data for performance
	}
	// Delete additionalMu
	if (additionalMu){ mxFree(additionalMu); }
	return 1;
}



/* The gateway function */
/****************************
%
%
% Matlab call:
%	[xytFPMat] = drugowitsch2Race(subInd, delta,  meshSize, mu, k, totalTime, dt, s0)
%
% Description:
%	Computes a 3D propagation of a diffusion particle into a 2D format (Please reshape outside of C since it's 
%	a lot easier to reshape a matrix in MATLAB). Using Method of Images (MoI) from "Family of closed-form solutions
%	for two-dimensional correlated diffusion processes" by Haozhe Shan, Ruben Moreno-Bote, and Jan Drugowitsch.
%	*************************************** STILL IN TESTING PHASE **********************************************
%	
%	Inputs:		subInd -> a 2xN vector with meshgrid points for propagation
%				delta -> The amount of change between mesh points in each dimension (X & Y), it should be identical for both
%				meshSize -> length(yi)*length(xi)
%				mu -> drift rate	
%				s0 -> Starting point of S
%				ti -> total Time
%				deltaT -> delta time (or time steps)
%				k -> parameter for correlation and number of images needed
%
%
*****************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Definition of mexFunction parameters
	nlhs = number of output arguments
	plhs = Array of output arguments
	nrhs = Number of input arguments
	prhs = Array of input arguments */

	/* variable declarations here*/
	/* Inputs */
	double* mu;
	double totalTime; double dt; int intTotalTime;
	double* subInd;
	double k;
	double* s0;
	double delta; 
	int xmeshSize;
	/*Temp variable values if not using Urgency Signal*/
	double urgencyMax = 0;
	double urgencyTauHalf = 0;
	/* Output(s) */
	//double* xytFPFinal = (double *)mxCalloc(xmeshSize*4000, sizeof(double));
	/* miscellenous */
	double err;

	/* code here */
	// First just make sure that you have the right number of inputs
	if (nrhs < 8)
		mexErrMsgTxt("Wrong number of arguments buddy \n"\
		" [xytFPMat] = drugowitsch2Race(subInd, delta, mu, k, totalTime, dt, s0)");
	/* Allocate memory for input parameters*/

	/* Simply load input parameters into C space*/
	// Get the drift rate & covariance matrix
	mu = mxGetPr(prhs[3]);
	// get total time
	totalTime = mxGetScalar(prhs[5]); 
	dt = mxGetScalar(prhs[6]);
	intTotalTime = (int)round(totalTime / dt + 1);
	// get meshSize
	xmeshSize = mxGetScalar(prhs[2]);

	/* Output(s) */
	double* xytFPFinal = (double *)mxCalloc(xmeshSize * intTotalTime, sizeof(double));

	//get S starting point
	s0 = mxGetPr(prhs[7]);
	//Get K, parameter for correlation and number of images needed
	k = mxGetScalar(prhs[4]);
	// Get Coordinates of MeshGrid (should be a Nx2 or 2xN)
	subInd = mxGetPr(prhs[0]);
	// Delta for Mesh Grid
	delta = mxGetScalar(prhs[1]);
	// Extra variables 
	double rho = -cos(pi / k);
	double alpha = ((k - 1) / k) * pi;
	double** globalSigma = (double **)mxCalloc(2, sizeof(double));

	for (int i = 0; i<2; ++i){
		globalSigma[i] = (double *)mxCalloc(2, sizeof(double));
		for (int j = 0; j<2; ++j){
			if (i == j) {
				globalSigma[i][j] = 1;
			}
			else {
				globalSigma[i][j] = rho;
			}
		}
	}
	// Load the urgency parameters if using them
	if (nrhs == 9){
		urgencyMax = mxGetScalar(prhs[8]);
	}
	if (nrhs == 10){
		urgencyMax = mxGetScalar(prhs[8]);
		urgencyTauHalf = mxGetScalar(prhs[9]);
	}
    
	/* Kiani at this point checks for memory allocation erros 
	but i think for now i will skip this, since I am not entirely 
	sure how he's doing this */

	/* Run program similar to kiani method in order to return error*/
	err = !drugo(subInd, delta, xmeshSize, mu, k, totalTime, dt, s0, alpha, globalSigma, urgencyMax, urgencyTauHalf, xytFPFinal);

	/* Prepare output data (which apparently is more difficult than one would think*/
	plhs[0] = mxCreateDoubleMatrix(xmeshSize, intTotalTime, mxREAL);
	mxFree(mxGetPr(plhs[0]));
	mxSetPr(plhs[0], xytFPFinal);

	/*return error message if error exist*/
	if (err){
		mexErrMsgTxt("Something went wrong with the MAIN code");
	}

}
