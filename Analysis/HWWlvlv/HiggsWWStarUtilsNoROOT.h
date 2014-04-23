#ifndef LESTER_HIGGSWWSTARUTILSNOROOT_H
#define LESTER_HIGGSWWSTARUTILSNOROOT_H

#include<cmath>
#include <cassert>

namespace Lester {

const bool DEBUG2 = false;

/// This is an implementation of what is written up in Lab Book 8B.  Roughly page 20. -- stop press -- no it isn't -- this is an implementation of a version that relaxes the mass constraint on side 2.
/// Higgs -> W (1) WStar (2)
/// W (1) -> s(vis) + p(invis) // here we do apply a W mass constraint
/// WStar (2) -> t(vis) + q(invis) // no W mass contraint applied
double HiggsWWStarMassLesterAtFixedKxKy_relaxedOnSide2Only(
	const double kx, const double ky,
	const double se, const double sx, const double sy, const double sz,
	const double te, const double tx, const double ty, const double tz,
        const double pmissx, const double pmissy,
	const double mw, 
	bool & wasSillyRef) {

wasSillyRef = true;

const double mwsq = mw*mw;
const double mssq = se*se-sx*sx-sy*sy-sz*sz;
// const double mtsq = te*te-tx*tx-ty*ty-tz*tz;

const double pkx1 = pmissx + kx;
const double pky1 = pmissy + ky;
// const double pkx2 = pmissx - kx;
// const double pky2 = pmissy - ky;

const double star1 = mwsq - mssq + sx*pkx1+sy*pky1;
// const double star2 = mwsq - mtsq + tx*pkx2+ty*pky2;

const double pk1sq = pkx1*pkx1 + pky1*pky1;
// const double pk2sq = pkx2*pkx2 + pky2*pky2;


const double etssq = se*se-sz*sz;
// const double ettsq = te*te-tz*tz;

bool good = true;

const double disc1 = star1*star1 - etssq*pk1sq;
// const double disc2 = star2*star2 - ettsq*pk2sq;

if (disc1<0) {
	good = false;
}

// where is invis 1?
const double px=0.5*(pmissx + kx);
const double py=0.5*(pmissy + ky);
const double pz_cat = (star1*sz + se*sqrt(disc1))/(2*etssq);
const double pz_dog = (star1*sz - se*sqrt(disc1))/(2*etssq);
const double pesq_cat = px*px + py*py + pz_cat*pz_cat;
const double pesq_dog = px*px + py*py + pz_dog*pz_dog;
const double pe_cat = sqrt(pesq_cat);
const double pe_dog = sqrt(pesq_dog);
 

// find 4-mom of everything else other than q:

const double alle_cat = pe_cat + se + te;
const double alle_dog = pe_dog + se + te;
const double allx     = px     + sx + tx;
const double ally     = py     + sy + ty;
const double allz_cat = pz_cat + sz + tz;
const double allz_dog = pz_dog + sz + tz;

// find mass of everything else other than q:

const double allmsq_cat = alle_cat*alle_cat - allx*allx - ally*ally - allz_cat*allz_cat;
const double allmsq_dog = alle_dog*alle_dog - allx*allx - ally*ally - allz_dog*allz_dog;

// find et of everything else other than q:

const double alletsq_cat = fabs(alle_cat*alle_cat - allz_cat*allz_cat); // have put the protection in, since alletsq_cat SHOULD be positive if the vectors have any pt.
const double alletsq_dog = fabs(alle_dog*alle_dog - allz_dog*allz_dog); // have put the protection in, since alletsq_dog SHOULD be positive if the vectors have any pt.
const double allet_cat = sqrt(alletsq_cat);
const double allet_dog = sqrt(alletsq_dog);

// find q_T
const double qx=0.5*(pmissx - kx);
const double qy=0.5*(pmissy - ky);
const double etq = sqrt(qx*qx + qy*qy);


// we know mass of q_T is zero by construction, so can now calculate the transverse mass of the all + q system treating "all" as visible (even though it contains a bit of p).
const double transversemasssq_cat = allmsq_cat + 0 + 2.0*( allet_cat*etq - allx*qx - ally*qy );
const double transversemasssq_dog = allmsq_dog + 0 + 2.0*( allet_dog*etq - allx*qx - ally*qy );

const double mhAAAsq = transversemasssq_cat;
const double mhBBBsq = transversemasssq_dog;

if (!good) {
        double power=1; // power to raise the -ve discriminant by in attemt to get out of bad region
	return pow(-disc1,power) + 14000;
}

const double mhsq = ((mhAAAsq <= mhBBBsq) ? mhAAAsq : mhBBBsq );

if (mhsq<0) {
	return -sqrt(-mhsq);
}

wasSillyRef = false;

return sqrt(mhsq);

} // end of function

/// This is an implementation of what is written up in Lab Book 8B.  Roughly page 20. -- stop press -- no it isn't -- this is an implementation of a version that relaxes the mass constraint on side 2.
// BOTH
/// Higgs -> W (1) WStar (2)
/// W (1) -> s(vis) + p(invis) // here we do apply a W mass constraint
/// WStar (2) -> t(vis) + q(invis) // no W mass contraint applied
// AND
/// Higgs -> W (1) WStar (2)
/// W (1) -> s(vis) + p(invis) // no W mass constraint applied
/// WStar (2) -> t(vis) + q(invis) // here we do apply a W mass contraint
// are considered, and the smaller is returned.

double HiggsWWStarMassLesterAtFixedKxKy(
	const double kx, const double ky,
	const double se, const double sx, const double sy, const double sz,
	const double te, const double tx, const double ty, const double tz,
        const double pmissx, const double pmissy,
	const double mw, 
	bool & wasSillyRef) {

  wasSillyRef = true;

  bool wasSilly1;
  bool wasSilly2;

  const double mh1 = HiggsWWStarMassLesterAtFixedKxKy_relaxedOnSide2Only(+kx,+ky,
									 se,sx,sy,sz,
									 te,tx,ty,tz,
									 pmissx,pmissy,
									 mw,
									 wasSilly1);
  const double mh2 = HiggsWWStarMassLesterAtFixedKxKy_relaxedOnSide2Only(-kx,-ky,
									 te,tx,ty,tz,
									 se,sx,sy,sz,
									 pmissx,pmissy,
									 mw,
									 wasSilly2);

  if (wasSilly1 && wasSilly2) {
    return (mh1<=mh2) ? mh1 : mh2;
  }
  if (wasSilly1 && !wasSilly2) {
    wasSillyRef = false;
    return mh2;
  }
  if (wasSilly2 && !wasSilly1) {
    wasSillyRef = false;
    return mh1;
  }
  if (!wasSilly2 && !wasSilly1) {
    wasSillyRef = false;
    return (mh1<=mh2) ? mh1 : mh2;
  }
  return -20; // error code
}

} // end of namespace Lester

#endif
