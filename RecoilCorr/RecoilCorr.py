import os
import ROOT
import math
import json
import glob
from array import array

inputsdir = {"2022preEE_NLO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_NLO/",
             "2022preEE_LO":    "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_LO/",
             "2022preEE_NNLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_NNLO/",
             "2022postEE_NLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_NLO/",
             "2022postEE_LO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_LO/",
             "2022postEE_NNLO": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_NNLO/",
             "2023preBPix_NLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_NLO/",
             "2023preBPix_LO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_LO/",
             "2023preBPix_NNLO": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_NNLO/",
             "2023postBPix_NLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_NLO/",
             "2023postBPix_LO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_LO/",
             "2023postBPix_NNLO": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_NNLO/",

             "2022preEE_NLO_appl":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_NLO_appl/",
             "2022preEE_LO_appl":    "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_LO_appl/",
             "2022preEE_NNLO_appl":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_RecoilCorrV2_NNLO_appl/",
             "2022postEE_NLO_appl":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_NLO_appl/",
             "2022postEE_LO_appl":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_LO_appl/",
             "2022postEE_NNLO_appl": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_RecoilCorrV2_NNLO_appl/",
             "2023preBPix_NLO_appl":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_NLO_appl/",
             "2023preBPix_LO_appl":    "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_LO_appl/",
             "2023preBPix_NNLO_appl":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_RecoilCorrV2_NNLO_appl/",
             "2023postBPix_NLO_appl":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_NLO_appl/",
             "2023postBPix_LO_appl":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_LO_appl/",
             "2023postBPix_NNLO_appl": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_RecoilCorrV2_NNLO_appl/",
            }
categories = ["0j_ZpT10",  "1j_ZpT10",  "2j_ZpT10",
              "0j_ZpT20",  "1j_ZpT20",  "2j_ZpT20",
              "0j_ZpT30",  "1j_ZpT30",  "2j_ZpT30",
              "0j_ZpT50",  "1j_ZpT50",  "2j_ZpT50",
              "0j_ZpTInf", "1j_ZpTInf", "2j_ZpTInf",
              "0j",        "1j",        "2j"]
finalstates = ["mumu", "ee"]
FinelyBinned = False

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

# Get all files containing histograms, and list of processes
hists = {}
fits = {}
Results = {}
comps = ["upara", "uperp", "response", "resolution"]
inputs = {}
for fs in finalstates:
  Results[fs] = {}
  inputs[fs] = {}
  thisdata = "Muon" if fs=="mumu" else "EGamma"
  for era in inputsdir:
    inputs[fs][era] = []
    for d in os.walk(inputsdir[era]):
      inputs[fs][era] += [os.path.join(d[0],f).replace(inputsdir[era],"") for f in d[2] if f.endswith(".root") and f"CR-{fs}" in f]
  processes = []
  for i in [f for era in inputs[fs] for f in inputs[fs][era] if f.endswith(f"_{categories[-1]}.root")]:
    i = i.split("/")[-1]
    cut = ["_HTauTau_"+era.split("_")[0] for era in inputsdir]
    cut = list(dict.fromkeys([c for c in cut if c in i]))
    assert len(cut)==1, f"Cuts are: {cut}"
    proc = i.split(cut[0])[0]
    if proc not in processes: processes.append(proc)
  
  def GetFitResult(fit):
    string = str(fit.GetExpFormula())
    for p in range(fit.GetNpar()):
      pname = fit.GetParName(p)
      string = string.replace(f"[{pname}]", str(fit.GetParameter(p)))
    string = string.replace("+-", "-")
    string = string.replace("--", "+")
    return string

  hists[fs] = {}
  fits[fs] = {}
  for era2 in inputsdir:
    era = era2[:-5] if era2.endswith("appl") else era2
    CorrectionApplied = era2.endswith("appl")
    hists[fs][era2] = {}
    for c in categories:
      if fs=="ee" and "ZpT" in c: continue
      hists[fs][era2][c] = {}
      for f in inputs[fs][era2]:
        if not f.endswith(f"_{c}.root"): continue
        proc = ""
        for p in processes:
          if f.split("/")[-1].startswith(p):
            proc = p
            break
        rfile = ROOT.TFile(os.path.join(inputsdir[era2], f), "read")
        if proc.startswith(thisdata): proc = thisdata # Combine data eras
        if proc.startswith("DY") and "10to50" not in proc: proc = "DY" # Combine DY, but treat 10-50 sample as background
        for comp in comps:
          key = [key.GetName() for key in rfile.GetListOfKeys() if key.GetName().startswith(comp) and ("fine" in key.GetName())==FinelyBinned]
          assert len(key)==1, f"Len is {len(key)}, Keys are {key}"
          h = rfile.Get(key[0])
          h.SetDirectory(0)
          if proc not in hists[fs][era2][c]:
            hists[fs][era2][c][proc]={}
          if comp not in hists[fs][era2][c][proc]:
            hists[fs][era2][c][proc][comp] = h
            hists[fs][era2][c][proc][comp].SetDirectory(0)
          else:
            hists[fs][era2][c][proc][comp].Add(h)
        rfile.Close()
  
    # Make DataMinusMC
    for c in hists[fs][era2]:
      hists[fs][era2][c]["DataMinusMC"] = {}
      for comp in comps:
        hists[fs][era2][c]["DataMinusMC"][comp] = hists[fs][era2][c][thisdata][comp].Clone()
      for proc in hists[fs][era2][c]:
        if proc in ["DY", thisdata, "DataMinusMC"]: continue
        for comp in hists[fs][era2][c][proc]:
          hists[fs][era2][c]["DataMinusMC"][comp].Add(hists[fs][era2][c][proc][comp], -1)

    # Rebin
    def RebinHist(histdy, histdata):
      binning = []
      for i in range(histdy.GetNbinsX()):
        binning.append(histdy.GetXaxis().GetBinLowEdge(i+1))
      binning.append(histdy.GetXaxis().GetBinUpEdge(i+1))
      for i in range(int(histdy.GetNbinsX()/2)):
        if histdy.GetBinContent(i+1) <= 0 or histdata.GetBinContent(i+1) <= 0:
          del binning[i+1]
          newhistdy = histdy.Rebin(len(binning)-1, histdy.GetName(), array("d", binning))
          newhistdata = histdata.Rebin(len(binning)-1, histdata.GetName(), array("d", binning))
          return RebinHist(newhistdy,newhistdata)
      for i in range(histdy.GetNbinsX(), int(histdy.GetNbinsX()/2), -1):
        if histdy.GetBinContent(i) <= 0 or histdata.GetBinContent(i) <= 0:
          del binning[i-1]
          newhistdy = histdy.Rebin(len(binning)-1, histdy.GetName(), array("d", binning))
          newhistdata = histdata.Rebin(len(binning)-1, histdata.GetName(), array("d", binning))
          return RebinHist(newhistdy,newhistdata)
      return histdy, histdata
    for c in hists[fs][era2]:
      for comp in comps:
        hists[fs][era2][c]["DY"][comp],hists[fs][era2][c]["DataMinusMC"][comp] = RebinHist(hists[fs][era2][c]["DY"][comp],hists[fs][era2][c]["DataMinusMC"][comp])
          

    # Remove negative bins
    for c in hists[fs][era2]:
      for proc in ["DataMinusMC", "DY"]:
        for comp in comps:
          for i in range(hists[fs][era2][c][proc][comp].GetNbinsX()):
            if hists[fs][era2][c][proc][comp].GetBinContent(i+1) < 0:
              hists[fs][era2][c][proc][comp].SetBinContent(i+1, 0.0)
  
    # Make fits
    for c in hists[fs][era2]:
      for comp in comps:
        for proc in ["DataMinusMC", "DY"]:
          hists[fs][era2][c][proc][comp].GetYaxis().SetTitle("Events")
          title = f"DY {era.split('_')[-1]} MC" if proc=="DY" else "Data minus non-DY MC"
          hists[fs][era2][c][proc][comp].SetName("Hist_"+era+"_"+c+"_"+comp+"_"+proc)
          hists[fs][era2][c][proc][comp].SetTitle(title)
        if fs=="mumu" and ("res" not in comp) and ("ZpT" in c) and not CorrectionApplied: # Only fit for upara/uperp
          if era not in Results[fs]: Results[fs][era] = {}
          if era not in fits[fs]: fits[fs][era] = {}
          if c not in Results[fs][era]: Results[fs][era][c] = {}
          if c not in fits[fs][era]:fits[fs][era][c] = {}
          Results[fs][era][c][comp] = {}
          fits[fs][era][c][comp] = {}
          for proc in ["DataMinusMC", "DY"]:
            Results[fs][era][c][comp][proc] = {}
            fitname = f"prefit_{era}_{c}_{proc}_{comp}"
            prefit = ROOT.TF1(fitname, "gaus", hists[fs][era2][c][proc][comp].GetXaxis().GetXmin(), hists[fs][era2][c][proc][comp].GetXaxis().GetXmax())
            hists[fs][era2][c][proc][comp].Fit(fitname)
            # Make double gaussian fit
            fitname = f"fit_{era}_{c}_{proc}_{comp}"
            if comp=="uperp":
              # UPerp: Two symmetrical gausians with mean at 0
              fit = ROOT.TF1(fitname, "[0]*exp(-0.5*(x/[1])^2)+[2]*exp(-0.5*(x/[3])^2)", hists[fs][era2][c][proc][comp].GetXaxis().GetXmin(), hists[fs][era2][c][proc][comp].GetXaxis().GetXmax())
              fit.SetParameter(0,prefit.GetParameter("Constant")*0.5)
              fit.SetParameter(2,prefit.GetParameter("Constant")*0.5)
              fit.SetParameter(1,prefit.GetParameter("Sigma")/math.sqrt(2))
              fit.SetParameter(3,prefit.GetParameter("Sigma")*math.sqrt(2))
            elif comp=="upara":
              # UPara: Two assymmetrical gausians
              fit = ROOT.TF1(fitname, "[0]*exp(-0.5*((x-[1])/[2])^2)*(x<[1])+[0]*exp(-0.5*((x-[1])/[3])^2)*(x>=[1])+[4]*exp(-0.5*((x-[5])/[6])^2)*(x<[5])+[4]*exp(-0.5*((x-[5])/[7])^2)*(x>=[5])", hists[fs][era2][c][proc][comp].GetXaxis().GetXmin(), hists[fs][era2][c][proc][comp].GetXaxis().GetXmax())
              fit.SetParameter(0,prefit.GetParameter("Constant")*0.5)
              fit.SetParameter(4,prefit.GetParameter("Constant")*0.5)
              fit.SetParLimits(0, 0.01*prefit.GetParameter("Constant"), 2*prefit.GetParameter("Constant"))
              fit.SetParLimits(4, 0.01*prefit.GetParameter("Constant"), 2*prefit.GetParameter("Constant"))
              fit.SetParameter(1,prefit.GetParameter("Mean"))
              fit.SetParameter(5,prefit.GetParameter("Mean"))
              fit.SetParameter(2,prefit.GetParameter("Sigma")/math.sqrt(2))
              fit.SetParameter(3,prefit.GetParameter("Sigma")/math.sqrt(2))
              fit.SetParameter(6,prefit.GetParameter("Sigma")*math.sqrt(2))
              fit.SetParameter(7,prefit.GetParameter("Sigma")*math.sqrt(2))
            status = hists[fs][era2][c][proc][comp].Fit(fitname, "S")
            if status.Status()!=0:
              print(">>>>> This fit did NOT converge!")
            # Normalize fit result to 1 before saving
            if comp=="uperp":
              norm = fit.GetParameter(0)*math.sqrt(2*math.pi)*abs(fit.GetParameter(1)) + fit.GetParameter(2)*math.sqrt(2*math.pi)*abs(fit.GetParameter(3))
              fit.SetParameter(0, fit.GetParameter(0)/norm)
              fit.SetParError(0, fit.GetParError(0)/norm)
              fit.SetParameter(2, fit.GetParameter(2)/norm)
              fit.SetParError(2, fit.GetParError(2)/norm)
            elif comp=="upara":
              norm = fit.GetParameter(0)*math.sqrt(2*math.pi)*(abs(fit.GetParameter(2))+abs(fit.GetParameter(3)))/2 + fit.GetParameter(4)*math.sqrt(2*math.pi)*(abs(fit.GetParameter(6))+abs(fit.GetParameter(7)))/2
              fit.SetParameter(0, fit.GetParameter(0)/norm)
              fit.SetParError(0, fit.GetParError(0)/norm)
              fit.SetParameter(4, fit.GetParameter(4)/norm)
              fit.SetParError(4, fit.GetParError(4)/norm)
            Results[fs][era][c][comp][proc]["Fit"] = GetFitResult(fit)
            Results[fs][era][c][comp][proc]["Histogram"] = [hists[fs][era2][c][proc][comp].GetBinContent(i) for i in range(hists[fs][era2][c][proc][comp].GetNbinsX()+2)]
            Results[fs][era][c][comp][proc]["Mean"] = prefit.GetParameter("Mean")
            Results[fs][era][c][comp][proc]["Sigma"] = prefit.GetParameter("Sigma")
            canv = ROOT.TCanvas(f'c_{era}_{c}_{proc}_{comp}', f'c_{era}_{c}_{proc}_{comp}', 600, 600)
            canv.cd()
            xname = "U_{parallel}" if comp=="upara" else "U_{perpendicular}"
            hists[fs][era2][c][proc][comp].GetXaxis().SetTitle(xname)
            hists[fs][era2][c][proc][comp].Draw("")
            canv.SetTopMargin(0.075)
            canv.SetBottomMargin(0.15)
            canv.SetLeftMargin(0.15)
            canv.SetRightMargin(0.05)
            fit.Draw("SAME")
            fits[fs][era][c][comp][proc] = fit
            if not os.path.exists(f"plots_{era}"): os.makedirs(f"plots_{era}")
            if not os.path.exists(f"plots_{era}/{c}"): os.makedirs(f"plots_{era}/{c}")
            canv.SaveAs(f"plots_{era}/{c}/RecoilCorrFit_{era}_{c}_{proc}_{comp}.png")
            binning = [hists[fs][era2][c][proc][comp].GetBinLowEdge(i+1) for i in range(hists[fs][era2][c][proc][comp].GetNbinsX())]
            binning.append(binning[-1] + hists[fs][era2][c][proc][comp].GetBinWidth(hists[fs][era2][c][proc][comp].GetNbinsX()))
            Results[fs][era][c][comp]["HistBinning"] = binning
        elif "res" in comp and CorrectionApplied:
          if era not in Results[fs]: Results[fs][era] = {}
          if c not in Results[fs][era]: Results[fs][era][c] = {}
          Results[fs][era][c][comp] = {}
          if comp=="response" and fs=="mumu" and ("ZpT" in c): # Get plain mean response from MC
            proc = "DY"
            Results[fs][era][c][comp][proc] = {"ResponseMean" : hists[fs][era2][c][proc][comp].GetMean()}
            print("Got Response Mean for:",fs,era,c,proc,comp,"=",hists[fs][era2][c][proc][comp].GetMean())
          if fs=="ee" and ("ZpT" not in c):
            if comp=="response":
              Results[fs][era][c][comp][proc] = {"ResponseUnc" : hists[fs][era2][c]["DataMinusMC"][comp].GetMean()-hists[fs][era2][c]["DY"][comp].GetMean()}
              print("Got Response Uncertainty for:",fs,era,c,proc,comp,"=",hists[fs][era2][c]["DataMinusMC"][comp].GetMean()-hists[fs][era2][c]["DY"][comp].GetMean())
            if comp=="resolution":
              Results[fs][era][c][comp][proc] = {"ResolutionUnc" : (hists[fs][era2][c]["DataMinusMC"][comp].GetRMS()-hists[fs][era2][c]["DY"][comp].GetRMS())/hists[fs][era2][c]["DY"][comp].GetRMS()}
              print("Got Resolution Uncertainty for:",fs,era,c,proc,comp,"=",(hists[fs][era2][c]["DataMinusMC"][comp].GetRMS()-hists[fs][era2][c]["DY"][comp].GetRMS())/hists[fs][era2][c]["DY"][comp].GetRMS())

print("=========================")
print("RESULTS")
print("=========================")
for fs in Results:
  for era in Results[fs]:
    for c in Results[fs][era]:
      for comp in Results[fs][era][c]:
        for proc in Results[fs][era][c][comp]:
          if proc=="HistBinning": continue
          print(f"{fs}, {era}, {c}, {proc}, {comp}:")
          for what in Results[fs][era][c][comp][proc]:
            if what=="Histogram": continue
            print(what,":",Results[fs][era][c][comp][proc][what])
          print("----------") 

with open('RecoilCorrFitResults.json', 'w') as file:
  json.dump(Results, file)

file = ROOT.TFile.Open("RecoilCorrFitResults.root", "RECREATE")
for fs in hists:
  for era in hists[fs]:
    for c in hists[fs][era]:
      for comp in comps:
        for proc in ["DataMinusMC", "DY"]:
          file.WriteObject(hists[fs][era][c][proc][comp], "Hist_"+fs+"_"+era+"_"+c+"_"+comp+"_"+proc)
          if "res" not in comp and "ZpT" in c and "appl" not in era:
            print("Writing fit for:",fs,era,c,comp,proc)
            file.WriteObject(fits[fs][era][c][comp][proc], "Fit_"+fs+"_"+era+"_"+c+"_"+comp+"_"+proc)
file.Close()
