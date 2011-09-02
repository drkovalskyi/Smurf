#ifndef LESTER_HIGGSWWSTARMASSBOUNDNOROOT_H
#define LESTER_HIGGSWWSTARMASSBOUNDNOROOT_H

#include <iostream>
#include "HiggsWWStarUtilsNoROOT.h"
#include "TRandom.h"
#include <cmath>

namespace Lester {

// A negative return value indicates an error!

double higgsWWStarMassBoundNoROOT(const double se, const double sx, const double sy, const double sz,
		      const double te, const double tx, const double ty, const double tz,
		      const double pmissx, const double pmissy,
		      const double mw) {

    double kxStart = 0;
    double kyStart = 0;
    double bestHMassSoFar = 0;

    bool haveValidStartPoint=false;
    bool bestPointWasSilly=false;

    // Attempt to get valid start point
    const double distFromWall /* AKA scan size */ = 2; // order of magnitude spread over which cauchy vals will be distributed.
   
    for (int i=0; i<10000; ++i) {

      // Warning! This doggedly tries to get a start point .. and may take forever when it is not possible start point

      const double theta = (gRandom->Uniform()-0.5)*TMath::Pi();
      const double distToStep = distFromWall*tan(theta); 
      const double angToStep = gRandom->Uniform()*TMath::Pi()*2.0;
      const double kx = distToStep * cos(angToStep);
      const double ky = distToStep * sin(angToStep);
      bool wasSilly;

      const double possHMass =
  	Lester::HiggsWWStarMassLesterAtFixedKxKy(kx,ky,
  						 se,sx,sy,sz,
  						 te,tx,ty,tz,
  						 pmissx,pmissy,
  						 mw,
  						 wasSilly);
      if (   (possHMass>=0 && !haveValidStartPoint) ||  (possHMass>=0 && haveValidStartPoint && possHMass < bestHMassSoFar)   )  {
  	bestHMassSoFar = possHMass;
  	haveValidStartPoint = true;
  	bestPointWasSilly = wasSilly;
  	kxStart = kx;
  	kyStart = ky;
  	if (Lester::DEBUG2) { 
  	    std::cout << "Found better start point with mass " << possHMass << std::endl;
  	}
  	if (!wasSilly) {
  	  // Don't need to work any harder ...
  	  break;
  	}
      }
    }
    if (Lester::DEBUG2) {
  	std::cout << "LESTER ANS " << bestHMassSoFar << "\tsilliness " << bestPointWasSilly << " valid " << haveValidStartPoint << std::endl;
    }

    if (haveValidStartPoint) {
      // Now we can attempt to minimise this function.
      
      double kxOld = kxStart;
      double kyOld = kyStart;

      double typicalStepSize = distFromWall;
      const double shrinkageFactor = 0.99;
      const double growthFactor = 2;
      
      while(typicalStepSize > 1.e-6) {
  	const double theta = (gRandom->Uniform()-0.5)*TMath::Pi();
  	const double distToStep = typicalStepSize*tan(theta); 
  	const double angToStep = gRandom->Uniform()*TMath::Pi()*2.0001;
  	const double newkx = kxOld + distToStep * cos(angToStep);
  	const double newky = kyOld + distToStep * sin(angToStep);
  	bool wasSilly;
   
  	const double possHMass =
  	  Lester::HiggsWWStarMassLesterAtFixedKxKy(newkx,newky,
  						   se,sx,sy,sz,
  						   te,tx,ty,tz,
  						   pmissx,pmissy,
  						   mw,
  						   wasSilly);
  	if ( possHMass>=0 && possHMass < bestHMassSoFar) {
  	  bestHMassSoFar = possHMass;
  	  bestPointWasSilly = wasSilly;
  	  kxOld = newkx;
  	  kyOld = newky;
  	  typicalStepSize *= growthFactor;
  	  if (Lester::DEBUG2) {
  		std::cout << "Found better point with mass " << possHMass << " while typical step " << typicalStepSize << std::endl;
  	  }
  	  
  	} // if point is an improvement
  	else {
  	  typicalStepSize *= shrinkageFactor;
  	} // point was not an improvement

      } // while we wantto keep going
      if (bestPointWasSilly) {
  	return -10;
      } else {
  	return bestHMassSoFar;
      }
    } // have valid start
    else {
      return -15;
    } // have invalid start
  }
} // end of namespace Lester
#endif
