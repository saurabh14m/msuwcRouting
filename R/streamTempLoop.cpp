#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix streamTempLoop(int timeLength, CharacterVector edgeIDs, IntegerVector orders, NumericMatrix velocities, NumericVector lengths, NumericMatrix RsurfSnow, NumericMatrix RsurfNoSnow, NumericMatrix flowqSub, NumericMatrix flowqOut, NumericMatrix flowqIn, NumericMatrix flowsRiv, NumericVector annualTmean, String by, List parentList, double K){

	NumericMatrix Twater(timeLength, edgeIDs.size())

	for(int timeStep=0; timeStep < timeLength; i++){

	  for(int i=0; i < edgeIDs.size(), i++){

		  int ord = orders[i];

		  double v = velocities(timeStep, i);

		  double len = lengths[i];

		  double qSnow = RsurfSnow(timeStep, i);
		  double qNoSnow = RsurfNoSnow(timeStep, i);
		  double qGw = flowqSub(timeStep, i);
		  double qIn = flowqIn(timeStep, i);
			  
		  double Tgw = annualTmean[timeStep] + 1.5;


		  if(timeStep > 0){
			  double TwaterOld = Twater(timeStep - 1, i);
			  double sRiv = flowsRiv(timeStep - 1, i);
		  } else {
			  double TwaterOld = 0;
			  double sRiv = 0;
		  }

		  double TairLocal = Tair(timeStep, i);


		  if(by == "month"){
			  double TairLag = TairLocal;
			  double lamda = 1.0;
		  } else {
			  //print(" by == day ??")
			  //TairLag //### Need to make some sort of lag for daily air temp
			  //lamda //###
		  }


		  double qLocal = (qSnow + qGw + qNoSnow);


		  double TwaterLocal = 0;

		  if(qLocal > 0){
			  TwaterLocal = ((Tsnow*qSnow)+(Tgw*qGw)+(lamda*TairLag)*qNoSnow)/qLocal;
		  }

		  double TwaterUpstream = 0;

		  if(ord > 1){
			
			  NumericVector parents = parentList[i];

			  double qInUpstream = 0;
			  double TWeight = 0;

			  if(parents.size() > 0){

				  for(int j=0; j < parents.size(); j++){

					  double qInUp = flowqOut(timeStep, parents[j]);

					  TWeight += Twater(timeStep, parents[j]) * qInUp;

					  qInUpstream += qInUp;
				  }

				  if(qInUpstream > 0){
					  TwaterUpstream = TWeight/qInUpstream;
				  }
			  }
		  }

		  double rSqSubTT;
		  double sRivTT;
		  double qInTT = 0;

		  if(len/v <= 1){

			rSqSubTT = len/(2*v);
			sRivTT = len/(2*v);

			if(ord > 1){
			  //##Get dimensionless factor, multiply times time to get TT?
			  qInTT  = len/v; //#gives fraction of timestep
			}

		  } else {
			//print("len/v > 1????")
			qInTT = 0;
			rSqSubTT = 1-v/(2*len);
			sRivTT = 1 - v/len;
		  }

		  double rSqSubCoeff = 1;
		  double sRivCoeff = 1;
		  double qInCoeff = 1;

		  if(rSqSubTT*K < 1){
			rSqSubCoeff = rSqSubTT*K;
		  }
		  if(sRivTT*K < 1){
			sRivCoeff = sRivTT*K;
		  }
		  if(qInTT*K < 1){
			qInCoeff = qInTT*K;
		  }



		  //reachTravelTime = simFlow$TT[timeStep, hydroID] ## TT (hour) travel time of water through the subbasain

		  //# Need to modify k ##

		  //TairLocal = 5

		  //#A function of travel time and shading, given a large enough travel time, it should be close to one

		  //if(TairLocal > 0){
		  if(true){
			  double TwaterLocalWarmed = TwaterLocal + (TairLocal - TwaterLocal)*rSqSubCoeff;
			  double TwaterQin = TwaterUpstream + (TairLocal - TwaterUpstream)*qInCoeff;
			  double TwaterSriv = TwaterOld + (TairLocal - TwaterOld)*sRivCoeff;

			//###OLD Twater code
			//Twater = TwaterInitial + (TairLocal - TwaterInitial) * K * reachTravelTime
		  } else {
			//epsilon = 2.5 //###Set to 2.5*C before sensitivity testing

			//TwaterLocalWarmed = TwaterLocal + ((TairLocal + epsilon) - TwaterLocal)*K*rSqSubTT

			//TwaterQin = TwaterUpstream + ((TairLocal + epsilon) - TwaterUpstream)*K*qInTT

			//TwaterSriv = TwaterOld + ((TairLocal + epsilon) - TwaterOld)*K*sRivTT
			//##Twater = TwaterInitial + ((TairLocal + epsilon) - TwaterInitial) * K * (reachTravelTime)
		  }

		  if(qIn == 0 && qLocal == 0 && sRiv == 0){
			double TwaterEdge = 0;
		  } else {
			double TwaterEdge = (TwaterQin*qIn + TwaterLocalWarmed*qLocal + TwaterSriv*sRiv)/(qIn+qLocal+sRiv);
		  }

		  if(TwaterEdge < 0){
			TwaterEdge = .1;
		  }

		  // Store values in results list
		  Twater(timeStep, i) = Twater;




		  //if(Twater > 50) stop()

	  }
	  }
	}
	return Twater
}
