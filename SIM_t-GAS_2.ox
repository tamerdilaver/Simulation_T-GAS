/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all Student's t GAS model parameters ( nu, omega, A and B)
**
**  Date:
**    	16/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.F. Hoogerheide & S.J. Koopman
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
	avThetaStar[0][2] = log(vTheta[2]);
	avThetaStar[0][3] = log(vTheta[3]-4)-log(100-vTheta[3]);
	avThetaStar[0][4] = log(vTheta[4]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adLambda, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega, const adLambda, const adGamma, const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);
	adOmega[0] = exp(vTheta[2]);
	adLambda[0] = 4+(100-4)*exp(vTheta[3])/(1+exp(vTheta[3]));
	adGamma[0] = exp(vTheta[4]);
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
	decl dAlpha, dBeta, dOmega, dLambda, dGamma;
	fGetPars( &dAlpha,  &dBeta, &dOmega, &dLambda, &dGamma, vTheta);

	decl dS2 = dGamma;					 											//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	//Please make these more efficient!!! delete the log()
	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = -1/2*log(M_PI)-1/2*log(dLambda-2) -1/2*log(dS2) -log(gammafact(dLambda/2))+log(gammafact((dLambda+1)/2)) -(dLambda+1)/2*log(1+ s_vY[i]^2 / ((dLambda-2)*dS2));
			//GAS recursion
			dS2 = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr( s_vY[i])/((dLambda-2)*  dS2))^(-1)*sqr( s_vY[i]) -   dS2) + dBeta*dS2;
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
	avTheta[0][2] = exp(vThetaStar[2]);
	avTheta[0][3] = 4+(100-4)*exp(vThetaStar[3])/(1+exp(vThetaStar[3]));
	avTheta[0][4] = exp(vThetaStar[4]);
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

fEstimateGAS(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adLambda_hat, const adGamma_hat){

	//initialise parameter values
	decl vTheta = zeros(iPARS,1);
	vTheta = <0.1 ; 0.99 ; 0.05 ; 7 ; 0.1>;			// Alpha, Beta, Omega, Lambda Startingvalues
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
	adOmega_hat[0] = vTheta[2];
	adLambda_hat[0] = vTheta[3];
	adGamma_hat[0] = vTheta[4];

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

		decl dAlpha_hat, dBeta_hat, dOmega_hat, dLambda_hat, dGamma_hat, vSE;
		vSE = fEstimateGAS(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dLambda_hat, &dGamma_hat);	 //Omega and lambda also estimated

		mTemp[i][0] =  (dAlpha_hat 	- dALPHA)	/vSE[0];	
		mTemp[i][1]	=  (dBeta_hat  	- dBETA)	/vSE[1];
		mTemp[i][2]	=  (dOmega_hat	- dOMEGA)	/vSE[2];
		mTemp[i][3]	=  (dLambda_hat - dLAMBDA)	/vSE[3];
		mTemp[i][4]	=  (dGamma_hat 	- dGAMMA)	/vSE[4];
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

fMonteCarlo2(const amAlpha, const amBeta, const amOmega, const amLambda, const amGamma, const amAlpha2, const amBeta2, const amOmega2, const amLambda2, const amGamma2){
	decl mTemp, mTempAlpha, mTempBeta, mTempOmega, mTempLambda, mTempGamma;
	decl mTemp2, mTempAlpha2, mTempBeta2, mTempOmega2, mTempLambda2, mTempGamma2;

	decl iPunten;
	iPunten = floor((iSIZE-1000)/iSTEPS);
	
	mTempAlpha = mTempBeta = mTempOmega  = mTempLambda = mTempGamma = zeros(iPunten,3);
	mTempAlpha2 = mTempBeta2 = mTempOmega2  = mTempLambda2 = mTempGamma2 = zeros(iPunten,3);
	mTemp = mTemp2 = zeros(iB,iPARS);

	decl iSize = iSIZE;
	
	for(decl j = 0; j<floor((iSize-1000)/iSTEPS) ; j++){
		iSIZE = 1000+iSTEPS*j;
		for(decl i = 0; i<iB ; i++){
			decl vReturns;
			fSimGAS2(dALPHA, dBETA, dOMEGA, dLAMBDA, dGAMMA, &vReturns, i);
	
			decl dAlpha_hat, dBeta_hat, dOmega_hat, dLambda_hat, dGamma_hat, vSE;
			vSE = fEstimateGAS(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dLambda_hat, &dGamma_hat);	 //Omega and Lambda also estimated
	
			mTemp[i][0] =  sqrt(iSIZE)*(dAlpha_hat - dALPHA);			//SQRT(T)*(\hat_\alpha_T - \alpha_0) ~ N(0, \SIGMA)
			mTemp[i][1]	=  sqrt(iSIZE)*(dBeta_hat  - dBETA);
			mTemp[i][2]	=  sqrt(iSIZE)*(dOmega_hat - dOMEGA);
			mTemp[i][3]	=  sqrt(iSIZE)*(dLambda_hat- dLAMBDA);
			mTemp[i][4]	=  sqrt(iSIZE)*(dGamma_hat-dGAMMA);

			//2nd part
			mTemp2[i][0] 	=  (dAlpha_hat - dALPHA)/vSE[0];				//(\hat_\alpha_T - \alpha_0)/SE(\hat_\alpha) ~ N(0, 1)
			mTemp2[i][1]	=  (dBeta_hat  - dBETA)/vSE[1];
			mTemp2[i][2]	=  (dOmega_hat - dOMEGA)/vSE[2];
			mTemp2[i][3]	=  (dLambda_hat- dLAMBDA)/vSE[3];
			mTemp2[i][4]	=  (dGamma_hat - dGAMMA)/vSE[4];
			
		}
		// v0.025_quantile, vMean, v0.975_quantile;				We get 95%-intervals

		decl vMeanTemp, vQ0025Temp, vQ0975Temp;
		
		vMeanTemp = meanc(mTemp);
		for(decl i=0;i<sizerc(vecindex(isdotnan(meanc(mTemp))));i++){
			vMeanTemp[vecindex(isdotnan(meanc(mTemp)))[i]] = meanc(deleter(mTemp[:][vecindex(isdotnan(meanc(mTemp)))[i]]));
		}
		vQ0025Temp = quantilec(mTemp,0.025);
		vQ0975Temp = quantilec(mTemp,0.975);

		mTempAlpha[j][0] = vQ0025Temp'[0];
		mTempAlpha[j][1] = vMeanTemp'[0];
		mTempAlpha[j][2] = vQ0975Temp'[0];
	
		mTempBeta[j][0] = vQ0025Temp'[1];
		mTempBeta[j][1] = vMeanTemp'[1];
		mTempBeta[j][2] = vQ0975Temp'[1];

		mTempOmega[j][0] = vQ0025Temp'[2];
		mTempOmega[j][1] = vMeanTemp'[2];
		mTempOmega[j][2] = vQ0975Temp'[2];
	
		mTempLambda[j][0] = vQ0025Temp'[3];
		mTempLambda[j][1] = vMeanTemp'[3];
		mTempLambda[j][2] = vQ0975Temp'[3];

		mTempGamma[j][0] = vQ0025Temp'[4];
		mTempGamma[j][1] = vMeanTemp'[4];
		mTempGamma[j][2] = vQ0975Temp'[4];

		vMeanTemp = meanc(mTemp2);
		for(decl i=0;i<sizerc(vecindex(isdotnan(meanc(mTemp2))));i++){
			vMeanTemp[vecindex(isdotnan(meanc(mTemp2)))[i]] = meanc(deleter(mTemp2[:][vecindex(isdotnan(meanc(mTemp2)))[i]]));
		}
		vQ0025Temp = quantilec(mTemp2,0.025);
		vQ0975Temp = quantilec(mTemp2,0.975);

		mTempAlpha2[j][0] = vQ0025Temp'[0];
		mTempAlpha2[j][1] = vMeanTemp'[0];	  //deletec()
		mTempAlpha2[j][2] = vQ0975Temp'[0];
	
		mTempBeta2[j][0] = vQ0025Temp'[1];
		mTempBeta2[j][1] = vMeanTemp'[1];
		mTempBeta2[j][2] = vQ0975Temp'[1];

		mTempOmega2[j][0] = vQ0025Temp'[2];
		mTempOmega2[j][1] = vMeanTemp'[2];
		mTempOmega2[j][2] = vQ0975Temp'[2];
	
		mTempLambda2[j][0] = vQ0025Temp'[3];
		mTempLambda2[j][1] = vMeanTemp'[3];
		mTempLambda2[j][2] = vQ0975Temp'[3];

		mTempGamma2[j][0] = vQ0025Temp'[4];
		mTempGamma2[j][1] = vMeanTemp'[4];
		mTempGamma2[j][2] = vQ0975Temp'[4];

	}

	amAlpha[0] 	= mTempAlpha;
	amBeta[0] 	= mTempBeta;
	amOmega[0] 	= mTempOmega;
	amLambda[0]	= mTempLambda;
	amGamma[0]	= mTempGamma;

	amAlpha2[0] = mTempAlpha2;
	amBeta2[0] 	= mTempBeta2;
	amOmega2[0] = mTempOmega2;
	amLambda2[0]= mTempLambda2;
	amGamma2[0]= mTempGamma2;

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
	dBETA = 0.99;
	dOMEGA = 0.05;
	dLAMBDA = 7;
	dGAMMA = 0.1;
	iPARS = 5;


///*
//** ..................................................................................	
//**	 		ASYMPTOTIC NORMALITY
//**	Get distributions of alpha and beta (to check for asymptotic normality)
//**..................................................................................
//*/

	//SET # OF SIMULATIONS 
	iB = 5000;			//5000 possible 
	iSIZE = 5000;		//5000 possible
	iSIMS = iB*iSIZE;
	vSTD_NORM = rann(iSIMS,1);
	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1, dLAMBDA);
	bSTD_ERROR = TRUE;

	decl dBeta_0;
	fGetBeta(iSIMS,  &dBeta_0);
	dBETA = dBeta_0;

	//DO MANY SIMULATIONS AND ESITMATIONS	
	decl mMonteCarlo;
	fMonteCarlo(&mMonteCarlo);	  

	//DRAW GRAPHS
	SetDrawWindow("SIM_t-GAS_2");

	DrawDensity(0, (mMonteCarlo[][0])', {"(i) Density $\hat A_1$ AsymN"});
	DrawDensity(1, (mMonteCarlo[][1])', {"(ii) Density $\hat B_1$ AsymN"});
	DrawDensity(2, (mMonteCarlo[][2])', {"(iii) Density $\hat\omega AsymN$"});
	DrawDensity(3, (mMonteCarlo[][3])', {"(iv) Density $\hat\lambda AsymN$"});
	DrawDensity(4, (mMonteCarlo[][4])', {"(v) Density $\hat \gamma AsymN$"});
	DrawDensity(0, vSTD_NORM', {"Z~N(0,1)"});
	DrawDensity(1, vSTD_NORM', {"Z~N(0,1)"});
	DrawDensity(2, vSTD_NORM', {"Z~N(0,1)"});
	DrawDensity(3, vSTD_NORM', {"Z~N(0,1)"});
	DrawDensity(4, vSTD_NORM', {"Z~N(0,1)"});
	
	DrawTitle(0,"(i) $\hat A_1$ AsymN");	
	DrawTitle(1,"(ii) $\hat B_1$ AsymN");
	DrawTitle(2,"(iii) $\hat \omega$ AsymN");
	DrawTitle(3,"(iv) $\hat \lambda$ AsymN");
	DrawTitle(4,"(v) $\hat \gamma$ AsymN");


	print("\nFirst Graph Finished at ",time(),"\n");
/*
** ..................................................................................	
**	 			CONSISTENCY
**	Check consistency for alpha and beta
** ..................................................................................
*/	

	//SET # OF SIMULATIONS 
	iB = 100;			  //100
	iSIZE = 100000;		  //10000
	iSIMS = iB*iSIZE;
	mSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iB,iSIZE,dLAMBDA);
	vSTD_STUDENT_T = vec(mSTD_STUDENT_T);	 //draw standardized student's t distributed 
	bSTD_ERROR = TRUE;

	//GET BETA SUCH THAT WE GET EQUALITY 
	//decl dBeta_0;
	fGetBeta(iSIMS,  &dBeta_0);
	dBETA = dBeta_0;
		
	//DO MANY SIMULATIONS AND ESITMATIONS
	decl mAlpha, mBeta, mOmega, mLambda, mGamma, mAlpha2, mBeta2, mOmega2, mLambda2, mGamma2;
	iSTEPS = iSIZE/20;				 	//steps of iSIZE/100 takes a while (steps of iSIZE/10 is faster)
	fMonteCarlo2(&mAlpha, &mBeta, &mOmega, &mLambda, &mGamma, &mAlpha2, &mBeta2, &mOmega2, &mLambda2, &mGamma2);

	//DRAW GRAPH

	Draw(5, mAlpha',1000,iSTEPS);
	Draw(6, mBeta',1000,iSTEPS);
	Draw(7, mOmega',1000,iSTEPS);
	Draw(8, mLambda',1000,iSTEPS);
	Draw(9, mGamma',1000,iSTEPS);
	DrawTitle(5,"(vi) A Cons");	
	DrawTitle(6,"(vii) B Cons");
	DrawTitle(7,"(viii) omega Cons");	
	DrawTitle(8,"(ix) lambda Cons");
	DrawTitle(9,"(x) gamma Cons");

	print("\nSecond Graph Finished at ",time(),"\n");


	Draw(10, mAlpha2',1000,iSTEPS);
	Draw(11, mBeta2',1000,iSTEPS);
	Draw(12, mOmega2',1000,iSTEPS);
	Draw(13, mLambda2',1000,iSTEPS);
	Draw(14, mGamma2',1000,iSTEPS);
	DrawTitle(10,"(xi) $\hat A_1$ NormC");	
	DrawTitle(11,"(xii) $\hat B_1$ NormC");
	DrawTitle(12,"(xiii) $\hat \omega$ NormC");	
	DrawTitle(13,"(xiv) $\hat \lambda$ NormC");
	DrawTitle(14,"(xv) $\hat \gamma$ NormC");
	ShowDrawWindow();
	print("\nThird Graph Finished at ",time(),"\n");
}