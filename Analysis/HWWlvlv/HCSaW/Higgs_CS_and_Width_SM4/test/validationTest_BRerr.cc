#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "HiggsCSandWidthSM4.cc"

using namespace std;

int main()
{

  ofstream fileOut;
  char* fileName_[4] = {"Htautau.txt","Hgamgam.txt","HWW.txt","HZZ.txt"};

  string channel[4] = {"H->tautau","H->gamgam","H->WW","H->ZZ"};

  HiggsCSandWidthSM4 *myCSW = new HiggsCSandWidthSM4();

  int IDs[4] = {2,8,10,11};
  
  for( int i = 0; i < 4; i++)
    {
      
      fileOut.open(fileName_[i]);
      fileOut << " mH   BR("+channel[i]+")   BR_Hff    BR_HVV    BR_Hgg    BR_Hgamgam " << endl;

      double mH;
      double BRChan;
      double sqrts = 7;
      bool spline = true;
      double BR_Hff, BR_HVV, BR_Hgg, BR_Hgamgam;

      for( double j = 110; j < 140; j += 0.5)
	{

	  mH = j;
	  BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);

	  fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV 
		  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;

	}



      for( double k = 140; k < 160; k++ )
	{

          mH = k;
	  BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);


	  fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV 
		  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;


	}

      for( double l = 160; l < 290; l += 2)
	{
          mH = l;
          BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);


          fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV
                  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;



	}

      for( double m = 290; m < 350; m += 5)
	{

          mH = m;
          BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);


          fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV
                  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;




	}

      for( double n = 350; n < 400; n += 10)
	{

          mH = n;
	  BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);


          fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV
                  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;




	}

      for( double q = 400; q <= 1000; q += 20 )
	{

          mH = q;
	  BRChan = myCSW->HiggsBR(IDs[i],mH,spline);
	  BR_Hff = myCSW->HiggsBRErr_Hff(IDs[i],mH,sqrts);
	  BR_HVV = myCSW->HiggsBRErr_HVV(IDs[i],mH,sqrts);
          BR_Hgg = myCSW->HiggsBRErr_Hgluglu(IDs[i],mH,sqrts);
          BR_Hgamgam = myCSW->HiggsBRErr_Hgamgam(IDs[i],mH,sqrts);


          fileOut << setw(6) << mH << " " << setw(6) << BRChan << "   " << setw(6) << BR_Hff << "   " << setw(6) << BR_HVV
                  << "   " << setw(6) << BR_Hgg << "   " << setw(6) << BR_Hgamgam << endl;



	}

      fileOut.close();

    }




  return 0;

}
