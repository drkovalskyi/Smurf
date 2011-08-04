#! /usr/bin/env python

def ComputeTriggerEfficiency(eta, pt, histo): 

    import math
    
    nbins = histo.GetXaxis().GetNbins()
    if pt > (histo.GetXaxis().GetBinLowEdge(nbins) + histo.GetXaxis().GetBinWidth(nbins)):
        pt = histo.GetXaxis().GetBinLowEdge(nbins) + (histo.GetXaxis().GetBinWidth(nbins)/2.0)

    nbins = histo.GetYaxis().GetNbins()
    if math.fabs(eta) > (histo.GetYaxis().GetBinLowEdge(nbins) + histo.GetYaxis().GetBinWidth(nbins)):
        eta = histo.GetYaxis().GetBinLowEdge(nbins) + (histo.GetYaxis().GetBinWidth(nbins)/2.0)

    binX = histo.GetXaxis().FindBin(pt)
    binY = histo.GetYaxis().FindBin(math.fabs(eta))
    return histo.GetBinContent(binX, binY)


def triggerEfficiency(lid, eta, pt, histos) :

    if lid[0] == 11 and lid[1] == 11: 
        eff_sgl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['single_e'])
        eff_sgl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['single_e'])
        eff_dbl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['double_e'])
        eff_dbl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['double_e'])
    
    elif lid[0] == 13 and lid[1] == 13:
        eff_sgl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['single_m'])
        eff_sgl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['single_m'])
        eff_dbl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['double_m'])
        eff_dbl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['double_m'])
    
    elif lid[0] == 13 and lid[1] == 11:
        eff_sgl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['single_m'])
        eff_sgl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['single_e'])
        eff_dbl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['cross_m'])
        eff_dbl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['cross_e'])
    
    elif lid[0] == 11 and lid[1] == 13:
        eff_sgl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['single_e'])
        eff_sgl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['single_m'])
        eff_dbl_1 = ComputeTriggerEfficiency(eta[0], pt[0], histos['cross_e'])
        eff_dbl_2 = ComputeTriggerEfficiency(eta[1], pt[1], histos['cross_m'])
    

    return eff_dbl_1*eff_dbl_2 + eff_sgl_2*(1-eff_dbl_1) + eff_sgl_1*(1-eff_dbl_2)



def leptonEfficiency(lid, eta, pt, histos):

    import math
    
    mypt   = min(pt,49.999)
    myeta  = min(math.fabs(eta),2.4999)
    prob = 1.0
    if lid == 13:
        ptbin = histos['muon'].GetXaxis().FindBin(mypt)
        etabin = histos['muon'].GetYaxis().FindBin(myeta)	 
        prob = histos['muon'].GetBinContent(ptbin,etabin)
        
    elif lid == 11:
        ptbin = histos['electron'].GetXaxis().FindBin(mypt)
        etabin = histos['electron'].GetYaxis().FindBin(myeta)	 
        prob = histos['electron'].GetBinContent(ptbin,etabin)
        
    return prob



def fakeRate(eta, pt, histo) :

    import math

    mypt   = min(pt,34.999)
    myeta  = min(math.fabs(eta),2.4999)

    # Si
    #ptbin = histo.GetXaxis().FindBin(mypt)
    #etabin = histo.GetYaxis().FindBin(myeta)	 
    #prob = histo.GetBinContent(ptbin,etabin)

    # Dima
    ptbin = histo.GetYaxis().FindBin(mypt)
    etabin = histo.GetXaxis().FindBin(myeta)	 
    prob = histo.GetBinContent(etabin,ptbin)

    
    if prob==1: return 1
    return prob/(1-prob)
    







