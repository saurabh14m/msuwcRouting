#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List streamTempLoop(int timeLength, CharacterVector edgeIDs, IntegerVector orders, NumericMatrix velocities, NumericVector lengths, NumericMatrix RsurfSnow, NumericMatrix RsurfNoSnow, NumericMatrix flowqSub, NumericMatrix flowqOut, NumericMatrix flowqIn, NumericMatrix flowsRiv, NumericVector annualTmean, String by, List parentList, double K, NumericMatrix Tair){

	NumericMatrix Twater(timeLength, edgeIDs.size());
	NumericMatrix TwUpstream(timeLength, edgeIDs.size());
	NumericMatrix TwLocal(timeLength, edgeIDs.size());

	double stepsLooped = 0;

	for(int timeStep=0; timeStep < timeLength; timeStep++){

	  for(int i=0; i < edgeIDs.size(); i++){

		  int ord = orders[i];

		  double v = velocities(timeStep, i);

		  double len = lengths[i];

		  double qSnow = RsurfSnow(timeStep, i);
		  double qNoSnow = RsurfNoSnow(timeStep, i);
		  double qGw = flowqSub(timeStep, i);
		  double qIn = flowqIn(timeStep, i);
			  
		  double Tgw = annualTmean[timeStep] + 1.5;

		  double TwaterOld = 0;
		  double sRiv = 0;

		  if(timeStep > 0){
			  TwaterOld = Twater(timeStep - 1, i);
			  sRiv = flowsRiv(timeStep - 1, i);
		  }

		  double TairLocal = Tair(timeStep, i);

		  double TairLag;
		  double lamda;

		  if(by == "month"){
			  TairLag = TairLocal;
			  lamda = 1.0;
		  } else {
			  //print(" by == day ??")
			  //TairLag //### Need to make some sort of lag for daily air temp
			  //lamda //###
		  }


		  double qLocal = (qSnow + qGw + qNoSnow);


		  double TwaterLocal = 0;

		  if(qLocal > 0){
			  TwaterLocal = ((0.1*qSnow)+(Tgw*qGw)+(lamda*TairLag)*qNoSnow)/qLocal;
		  }

		  double TwaterUpstream = 0;

		  if(ord > 1){
			
			  NumericVector parents = parentList[i];

			  double qInUpstream = 0;
			  double TWeight = 0;

			  if(parents.size() > 0){

				  for(int j=0; j < parents.size(); j++){

					  double qInUp = flowqOut(timeStep, parents[j]-1);

					  TWeight += Twater(timeStep, parents[j]-1) * qInUp;

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
			rSqSubTT = 1-v/(2*len);
			sRivTT = 1 - v/len;
		  }

		  double rSqSubCoeff = 1.0;
		  double sRivCoeff = 1.0;
		  double qInCoeff = 1.0;

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
		  //


		  //if(true){
			  double TwaterLocalWarmed = TwaterLocal + (TairLocal - TwaterLocal)*rSqSubCoeff;
			  double TwaterQin = TwaterUpstream + (TairLocal - TwaterUpstream)*qInCoeff;
			  double TwaterSriv = TwaterOld + (TairLocal - TwaterOld)*sRivCoeff;

			//###OLD Twater code
			//Twater = TwaterInitial + (TairLocal - TwaterInitial) * K * reachTravelTime
		  //} else {
			//epsilon = 2.5 //###Set to 2.5*C before sensitivity testing

			//TwaterLocalWarmed = TwaterLocal + ((TairLocal + epsilon) - TwaterLocal)*K*rSqSubTT

			//TwaterQin = TwaterUpstream + ((TairLocal + epsilon) - TwaterUpstream)*K*qInTT

			//TwaterSriv = TwaterOld + ((TairLocal + epsilon) - TwaterOld)*K*sRivTT
			//##Twater = TwaterInitial + ((TairLocal + epsilon) - TwaterInitial) * K * (reachTravelTime)
		  //}

		  double TwaterEdge = 0;

		  if(qIn != 0 || qLocal != 0 || sRiv != 0){
			TwaterEdge = (TwaterQin*qIn + TwaterLocalWarmed*qLocal + TwaterSriv*sRiv)/(qIn+qLocal+sRiv);
		  }

		  if(TwaterEdge < 0){
			TwaterEdge = .1;
		  }

		  //Rcpp::Rcout << "TwaterEdge: " << TwaterEdge << std::endl;

		  // Store values in results list
		  Twater(timeStep, i) = TwaterEdge;
		  TwUpstream(timeStep, i) = TwaterUpstream;
		  TwLocal(timeStep, i) = TwaterLocal;




		  //if(Twater > 50) stop()

		  stepsLooped++;
	  }
	}
	Rcpp::Rcout << "Steps Looped: " << stepsLooped << std::endl;
	return List::create(Twater, TwUpstream, TwLocal);
}
