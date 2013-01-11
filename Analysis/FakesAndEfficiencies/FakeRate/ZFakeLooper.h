
#ifndef ZFAKELOOPER_H
#define ZFAKELOOPER_H

#include "Enums.h"

// C++ includes
#include <iostream>
#include <vector>

#include "../../../Core/SmurfTree.h"
#include "TChain.h"
#include "TH1F.h"

class ZFakeLooper {
    public:
        ZFakeLooper();
        ~ZFakeLooper() {};
        int Loop(bool isData, TChain* chain, const char* name, const std::vector<unsigned int>& ptThresholds);

        void setLumi(float lumi);

    private:

        float lumi_;

};

#endif

