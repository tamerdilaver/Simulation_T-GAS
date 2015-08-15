/*
**	Master Thesis Econometrics
**
**  Purpose:
**  	For some fixed parameter values simulate and estimate T-GAS(1,1) model parameters
**		with Maximum Likelikhood many times. s.t. stationarity conditions of Bougerol's Theorem
**
**  Date:
**    	15/08/2015
**
**  Author:
**	  	Tamer Dilaver
**
**	Supervisor:
**		Fransisco Blasques
**
*/
														
#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iB;	 					//Repeats
static decl iSIZE;					//Size of time series
static decl iSTEPS;					//#Steps to divide the size
static decl iSIMS;					//# of Zt ~ N(0,1)
static decl dLAMBDA;				//Degrees of freedom
static decl dALPHA;					//dALPHA is actually dA (please change this later)
static decl dBETA;
static decl dOMEGA;				 
static decl dGAMMA;
static decl dOMEGA_ERROR;
static decl dGAMMA_ERROR;
static decl iPARS;					//number of parameters
static decl vSTD_NORM;				// Zt ~ N(0,1)
static decl vSTD_STUDENT_T;			// Zt ~ TID(lambda)
static decl mSTD_STUDENT_T;
static decl s_vY; 					//Simulated returns
static decl	bSTD_ERROR;				//0 or 1 binary


/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/
fTransform2(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);

	return 1;
}

/*
**  Function:	Transform parameters back
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/
fTransformBack2(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adBeta, vTheta
**
**  Output: 	1 
*/
fGetPars2(const adBeta, const vTheta){

	adBeta[0] = exp(vTheta[0]);
	return 1;
}

/*
**  Function:	Calculate targetvalue -[ E log(alpha_0 z_t^2 + beta_0) ]^2 for given parameters
**
**  Input:		vTheta [parametervalues], adFunc [adres functionvalue], avScore [the score],  amHessian [hessianmatrix]
**
**  Output:		1
**
*/
fExpectation(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanc(log(fabs(dALPHA*((dLAMBDA+3)/dLAMBDA) *((dLAMBDA+1)/(dLAMBDA-2)*(1+(vSTD_STUDENT_T.^2)/(dLAMBDA-2)).^(-1).* vSTD_STUDENT_T.^2 -1) + dBeta))))^2;
	return 1;
}

/*
**  Function:	Get the value for beta subject to Elog(alpha_0 (v+3)/v [(v+1)/(v+2)(1+z_t^2/(v-2))^-1 z_t^2 - 1] + beta_0) = 0 for given parameters	and simulated Zt's
**
**  Input:		iSims, adBeta  
**
**  Output:		1
**
*/
fGetBeta(const iSims, const adBeta)
{
	decl vTheta, vThetaStart, vThetaStar, dFunc, iA;

	//initialise startparameter(s)
	vTheta = zeros(1,1);
	vTheta = <0.9>;			// dBeta
	vThetaStart = vTheta;

	//transform startparameter(s)
	fTransform2(&vThetaStar, vTheta);

	//maximise
	iA=MaxBFGS(fExpectation, &vThetaStar, &dFunc, 0, TRUE);
	
	//Transform thetasStar back
   	fTransformBack2(&vTheta, vThetaStar);

	print("\nStart & Optimal parameter(s) with A_0 fixed at ",dALPHA,",lambda fixed at ",dLAMBDA," and ",iSIMS," simulations such that we get I(1). \n",
          "%r", { "dBeta"},
          "%c", {"thetaStart","theta"}, vThetaStart~vTheta);
	adBeta[0] = vTheta[0];

	return 1;
}

/*
**  Function:	Simulate GAS returns for given parameters
**
**  Input:		dAlpha, dBeta, dOmega, dLambda, avReturns, iIteration [to get different Zt's]
**
**  Output:		1
**
*/

fSimGAS(const dAlpha, const dBeta, const dOmega, const dLambda, const dGamma, const avReturns, const iIteration){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma; //by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_STUDENT_T[(i + (iIteration*iSIZE))];
		vH[i+1] = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr(vTemp[i])/((dLambda-2)*vH[i]))^(-1)*sqr(vTemp[i])-vH[i]) + dBeta*vH[i];
	}					  

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
	return 1;
}

fSimGAS2(const dAlpha, const dBeta, const dOmega, const dLambda, const dGamma, const avReturns, const iIteration){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma; //by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*mSTD_STUDENT_T[iIteration][i];
		vH[i+1] = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr(vTemp[i])/((dLambda-2)*vH[i]))^(-1)*sqr(vTemp[i])-vH[i]) + dBeta*vH[i];
	}					  

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
	return 1;
}

/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);
	avThetaStar[0][1] = log(vTheta[1]);						
	avThetaStar[0][2] = log(vTheta[2]-4)-log(100-vTheta[2]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adLambda, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adLambda,  const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);								  
	adLambda[0] = 4+(100-4)*exp(vTheta[2])/(1+exp(vTheta[2]));
	return 1;
}

/*
**  Function:	Calculates average value loglikelihood for GAS given parameter values
**
**  Input: 		vTheta [parametervalues], adFunc [adres functievalue], avScore [the score], amHessian [hessianmatrix]
**
**  Output:		1
**
*/

fLogLike_GAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dAlpha, dBeta, dLambda;
	fGetPars( &dAlpha, &dBeta, &dLambda, vTheta);

	decl dS2 = dGAMMA_ERROR;					 											//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	//Please make these more efficient!!! delete the log()
	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = -1/2*log(M_PI)-1/2*log(dLambda-2) -1/2*log(dS2) -log(gammafact(dLambda/2))+log(gammafact((dLambda+1)/2)) -(dLambda+1)/2*log(1+ s_vY[i]^2 / ((dLambda-2)*dS2));
			//GAS recursion
			dS2 = dOMEGA_ERROR + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr( s_vY[i])/((dLambda-2)*  dS2))^(-1)*sqr( s_vY[i]) -   dS2) + dBeta*dS2;
	}	
							 
	adFunc[0] = sumc(vLogEta)/sizerc(s_vY); 									 	//Average
	return 1;
}

/*
**  Function:	Transform parameters back
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/

fTransformBack(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	avTheta[0][1] = exp(vThetaStar[1]);									
	avTheta[0][2] = 4+(100-4)*exp(vThetaStar[2])/(1+exp(vThetaStar[2]));
	//actually need to restrict dLambda_hat between (4,100)
	//otherwise there will be no convergence for small samples that occur Gaussian
	return 1;
}

/*
**  Function:	calculate standard errors
**
**  Input: 		vThetaStar
**
**  Output: 	vStdErrors
*/

fSigmaStdError(const vThetaStar){

 		decl iN, mHessian, mHess, mJacobian, vStdErrors, vP;

		iN 			= sizerc(s_vY);
		Num2Derivative(fLogLike_GAS, vThetaStar, &mHessian);
		//NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		//mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		mHessian 	= invertgen(-iN*mHessian);
		

		decl vDelta = ones(iPARS,1);		 //Jacobiaan Analytisch
		
		vDelta[0]	=	exp(vThetaStar[0]);
		vDelta[1]	=	exp(vThetaStar[1]);
		vDelta[2]	=	exp(vThetaStar[2]);
		vDelta[3]	=	(100-4)*exp(vThetaStar[3])/(1+exp(vThetaStar[3]))^2;
		vDelta[4]	=	exp(vThetaStar[4]);

		vStdErrors 	= vDelta.*sqrt(diagonal(mHessian)');
		
		return 	vStdErrors;
}

/*
**  Function:	Estimate GAS parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adLambda_hat (dBeta_0 not necessary)
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimateGAS(const vReturns, const adAlpha_hat, const adBeta_hat, const adLambda_hat){

	//initialise parameter values
	decl vTheta = zeros(iPARS,1);
	vTheta = <0.1 ; 0.99 ; 7>;			// Alpha, Beta, Omega, Lambda Startingvalues
	decl vThetaStart = vTheta;

	//globalize returns and vectorize true pars
	s_vY = vReturns;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_GAS, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return alpha, beta, omega and lambda
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];		  
	adLambda_hat[0] = vTheta[2];	  

	if(bSTD_ERROR){	//only do this for fMonteCarlo2
		decl vSigmaStdError = fSigmaStdError(vThetaStar);
		return vSigmaStdError;
	}else{
		return 1; //otherwise return 1 and end function
	}
}

/*
**  Function:	Simulates and Estimates GAS data and parameters many times
**				to illustrate Asymptotic normality
**
**  Input: 		amMonteCarlo [matrix of many estimated parameters],  dBeta_0;
**
**  Output: 	1
*/

fMonteCarlo(const amMonteCarlo){
	decl mTemp;
	mTemp = zeros(iB,iPARS);

	for(decl i = 0; i<iB ; i++){
		decl vReturns;
		fSimGAS(dALPHA, dBETA, dOMEGA, dLAMBDA, dGAMMA, &vReturns, i);

		bSTD_ERROR = FALSE;
		decl dAlpha_hat, dBeta_hat, dLambda_hat;
		fEstimateGAS(vReturns, &dAlpha_hat, &dBeta_hat, &dLambda_hat);	 //Omega and lambda also estimated

			mTemp[i][0] =  (dAlpha_hat - dALPHA);	
			mTemp[i][1]	=  (dBeta_hat  - dBETA);
			mTemp[i][2]	=  (dLambda_hat- dLAMBDA);
	}
	amMonteCarlo[0] = mTemp;
	return 1;
}

/*
**  Function:	Simulated and Estimates GAS data and parameters many times
**				to illustrate consistency it returns minimum, mean and maximum values for the estimated parameters
**
**  Input: 		amAlpha [matrix containing the min, max and mean of estimated alpha],
**				amBeta [matrix containing the min, max and mean of estimated beta], 
**				amOmega [matrix containing the min, max and mean of estimated omega],
**				amLambda [matrix containing the min, max and mean of estimated lambda],	dBETA
**
**  Output: 	1
*/

fMonteCarlo2(const amAlpha, const amBeta, const amLambda, const amAlpha2, const amBeta2, const amLambda2){
	decl mTemp, mTempAlpha, mTempBeta, mTempLambda;
	decl mTemp2, mTempAlpha2, mTempBeta2, mTempLambda2;
	mTempAlpha = mTempBeta  = mTempLambda = zeros((iSIZE/iSTEPS),3);
	mTempAlpha2 = mTempBeta2  = mTempLambda2 = zeros((iSIZE/iSTEPS),3);
	mTemp = mTemp2 = zeros(iB,iPARS);

	decl iSize = iSIZE;

	for(decl j = 0; j<(iSize/iSTEPS) ; j++){
		iSIZE = ((iSTEPS)*(j+1));
		for(decl i = 0; i<iB ; i++){
			decl vReturns;
			fSimGAS2(dALPHA, dBETA, dOMEGA, dLAMBDA, dGAMMA, &vReturns, i);
	
			decl dAlpha_hat, dBeta_hat, dLambda_hat, vSE;
			vSE = fEstimateGAS(vReturns, &dAlpha_hat, &dBeta_hat, &dLambda_hat);	 //Omega and Lambda also estimated
	
			mTemp[i][0] =  sqrt(iSIZE)*(dAlpha_hat - dALPHA);			//SQRT(T)*(\hat_\alpha_T - \alpha_0) ~ N(0, \SIGMA)
			mTemp[i][1]	=  sqrt(iSIZE)*(dBeta_hat  - dBETA);  
			mTemp[i][2]	=  sqrt(iSIZE)*(dLambda_hat- dLAMBDA);

			//2nd part
			mTemp2[i][0] 	=  (dAlpha_hat - dALPHA)/vSE[0];				//(\hat_\alpha_T - \alpha_0)/SE(\hat_\alpha) ~ N(0, 1)
			mTemp2[i][1]	=  (dBeta_hat  - dBETA)/vSE[1];	 
			mTemp2[i][2]	=  (dLambda_hat- dLAMBDA)/vSE[2];
			
		}
		// v0.025_quantile, vMean, v0.975_quantile;				We get 95%-intervals
		mTempAlpha[j][0] = quantilec(mTemp[][],0.025)'[0];
		mTempAlpha[j][1] = meanc(mTemp[][])'[0];
		mTempAlpha[j][2] = quantilec(mTemp[][],0.975)'[0];
	
		mTempBeta[j][0] = quantilec(mTemp[][],0.025)'[1];
		mTempBeta[j][1] = meanc(mTemp[][])'[1];
		mTempBeta[j][2] = quantilec(mTemp[][],0.975)'[1];
		
		mTempLambda[j][0] = quantilec(mTemp[][],0.025)'[2];
		mTempLambda[j][1] = meanc(mTemp[][])'[2];
		mTempLambda[j][2] = quantilec(mTemp[][],0.975)'[2];
															
		//2nd part: v0.025_quantile, v0.5_quantile, v0.975_quantile;	
		mTempAlpha2[j][0] = quantilec(mTemp2[][],0.025)'[0];
		mTempAlpha2[j][1] = meanc(mTemp2[][])'[0];
		mTempAlpha2[j][2] = quantilec(mTemp2[][],0.975)'[0];
	
		mTempBeta2[j][0] = quantilec(mTemp2[][],0.025)'[1];
		mTempBeta2[j][1] = meanc(mTemp2[][])'[1];
		mTempBeta2[j][2] = quantilec(mTemp2[][],0.975)'[1];
															
		mTempLambda2[j][0] = quantilec(mTemp2[][],0.025)'[2];
		mTempLambda2[j][1] = meanc(mTemp2[][])'[2];
		mTempLambda2[j][2] = quantilec(mTemp2[][],0.975)'[2];
															
	}

	amAlpha[0] 	= mTempAlpha;
	amBeta[0] 	= mTempBeta;	 
	amLambda[0]	= mTempLambda;	 

	amAlpha2[0] = mTempAlpha2;
	amBeta2[0] 	= mTempBeta2;	 
	amLambda2[0]= mTempLambda2;	 

	return 1;
}

/*
**				MAIN PROGRAM
**
**  Purpose:	For a given alpha, omega and lambda, get beta such that E log(alpha_0 z_t^2 + beta_0) = 0.
**				Simulate GAS returns for alpha, omega, lambda and the found beta many time.
**				Estimate GAS parameters alpha, beta, omega and lambda.
**
**  Input: 		dALPHA, dOMEGA, dLAMBDA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/

main(){

	//SET PARAMETERS
	dALPHA = 0.1;
//	dBETA = 0.99;
	dOMEGA = 0.05;
	dLAMBDA = 7;
	dGAMMA = 0.1;
	iPARS = 3;

//	//SET ERROR FOR ALPHA AND GAMMA
//	dOMEGA_ERROR = dOMEGA/5;		 //must be strict larger then -	dOmega (i.e.  > -dOmega*999/1000)
//	dGAMMA_ERROR = dGAMMA/5;		 //must be strict larger then -	dGAMMA (i.e.  > -dGamma*999/1000)
//	if(dOMEGA_ERROR!=0 && dGAMMA_ERROR!=0){
//		print("\ndOMEGA_ERROR is ",dOMEGA_ERROR, " and dGAMMA_ERROR is ", dGAMMA_ERROR,"\n");
//	}


/*
** ..................................................................................	
**	 		Check Bias
**..................................................................................
*/

	//SET # OF SIMULATIONS 
	iB = 100; 			//max 5000
	iSIZE = 1000;		//max 5000
	iSIMS = iB*iSIZE;
	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1, dLAMBDA);
	//bSTD_ERROR = TRUE;

	//GET BETA SUCH THAT WE GET EQUALITY "Elog(alpha_0 z_t^2 + beta_0) = 0".
	decl dBeta_0;
	fGetBeta(iSIMS,  &dBeta_0);
	dBETA = dBeta_0;

	
	decl dOmega_min, dOmega_max, dOmega_step, dGamma_min, dGamma_max, dGamma_step;
	dOmega_min = dOMEGA*0.1; //max(0.001,dOMEGA*0.1);
	dOmega_max = dOMEGA*2;
	dOmega_step = dOMEGA/10;

	dGamma_min = dGAMMA*0.1; //max(0.01,dGAMMA*0.1);
	//print(dGamma_min,"\n");
	dGamma_max = dGAMMA*2;
	//print(dGamma_max,"\n");
	dGamma_step = dGAMMA/10;
	//print(dGamma_step,"\n");
	
//	decl vReturns;
//	fSimGARCH(dALPHA, dBETA, dOMEGA, dGAMMA, &vReturns, 0);

	decl vX, vY, mZ, mZ_2, mZ_3;
	
	vX = zeros((dOmega_max-dOmega_min)/dOmega_step,1); //vX are sigma2_eta
	vY = zeros((dGamma_max-dGamma_min)/dGamma_step,1);

//	print(sizerc(vX),"\n");
//	print(sizerc(vY),"\n");
	
	mZ = mZ_2 = mZ_3 = zeros(sizerc(vX),sizerc(vY));

	
	for(decl i = 0; i<sizerc(vX);i++){
		vX[i]  = (i+1)*dOmega_step;
		
		for(decl w = 0; w<sizerc(vY) ; w++){
			vY[w] = (w+1)*dGamma_step;
			dOMEGA_ERROR = vX[i];
			dGAMMA_ERROR = vY[w];

			decl dAlpha_hat, dBeta_hat, vSE;

			//DO MANY SIMULATIONS AND ESITMATIONS	
			decl mMonteCarlo;
			fMonteCarlo(&mMonteCarlo);
			
			mZ[i][w] 	= meanc(mMonteCarlo[][0]);
			mZ_2[i][w] 	= meanc(mMonteCarlo[][1]);
			mZ_3[i][w] 	= meanc(mMonteCarlo[][2]);
			//print(i~w~mZ[i][w]~mZ_2[i][w]);
		}
		 
	}

	decl iMaxCorrect = maxc(maxc(-fabs(mZ))');
	decl iMaxCorrect2 = maxc(maxc(-fabs(mZ_2))');
	decl iMaxCorrect3 = maxc(maxc(-fabs(mZ_3))');
	


	print("\nMin estimate error for A: ", iMaxCorrect);
	//print("\nWith $\hat\omega$ = ",vX[maxcindex(maxc( -mZ' )')], " $\hat\gamma$ = ", vY[maxcindex(maxc(-mZ)')]);

	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "omega",  "gamma"},
          	"%c", {"thetaStart","theta"}, (dOMEGA|dGAMMA)~(vX[maxcindex(maxc( -fabs(mZ)' )')]|vY[maxcindex(maxc(-fabs(mZ))')]));

	print("\nMin estimate error for B: ", iMaxCorrect2);
	//print("\nWith $\hat\omega$ = ",vX[maxcindex(maxc( -mZ' )')], " $\hat\gamma$ = ", vY[maxcindex(maxc(-mZ)')]);

	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "omega",  "gamma"},
          	"%c", {"thetaStart","theta"}, (dOMEGA|dGAMMA)~(vX[maxcindex(maxc( -fabs(mZ_2)' )')]|vY[maxcindex(maxc(-fabs(mZ_2))')]));
			
	print("\nMin estimate error for Nu: ", iMaxCorrect3);
	//print("\nWith $\hat\omega$ = ",vX[maxcindex(maxc( -mZ' )')], " $\hat\gamma$ = ", vY[maxcindex(maxc(-mZ)')]);

	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "omega",  "gamma"},
          	"%c", {"thetaStart","theta"}, (dOMEGA|dGAMMA)~(vX[maxcindex(maxc( -fabs(mZ_3)' )')]|vY[maxcindex(maxc(-fabs(mZ_3))')]));

	SetDrawWindow("SIM_t-GAS_5");
	DrawXYZ( 0,vY, vX, mZ,0 , 	"$\gamma$", "$\omega$", "$E(\hat A-A_0)$");
	DrawXYZ( 1,vY, vX, mZ_2,0 , "$\gamma$", "$\omega$", "$E(\hat B-B_0)$");
	DrawXYZ( 2,vY, vX, mZ_3,0 , "$\gamma$", "$\omega$", "$E(\hat nu-nu_0)$");
	DrawTitle(0,"(i) $E(\hat A-A_0)$ for fixed $\omega$ and $\gamma$");
	DrawTitle(1,"(ii) $E(\hat B-B_0)$ for fixed $\omega$ and $\gamma$");
	DrawTitle(2,"(iii) $E(\hat \lambda-\lambda_0)$ for fixed $\omega$ and $\gamma$");
	ShowDrawWindow();


}
