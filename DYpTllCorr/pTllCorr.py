import os
import ROOT
import math
import json
from array import array
import numpy as np

inputsdir = {"2022preEE_NLO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_PTLLCorrV2_NLO/",
             "2022preEE_LO":    "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_PTLLCorrV2_LO/",
             "2022preEE_NNLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022preEE_PTLLCorrV2_NNLO/",
             "2022postEE_NLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_PTLLCorrV2_NLO/",
             "2022postEE_LO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_PTLLCorrV2_LO/",
             "2022postEE_NNLO": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2022postEE_PTLLCorrV2_NNLO/",
             "2023preBPix_NLO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_PTLLCorrV2_NLO/",
             "2023preBPix_LO":    "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_PTLLCorrV2_LO/",
             "2023preBPix_NNLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023preBPix_PTLLCorrV2_NNLO/",
             "2023postBPix_NLO":  "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_PTLLCorrV2_NLO/",
             "2023postBPix_LO":   "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_PTLLCorrV2_LO/",
             "2023postBPix_NNLO": "/afs/cern.ch/work/d/dmroy/HTauTau/CMSSW_13_0_10/src/TauFW/Plotter/jobs/HTauTau_2023postBPix_PTLLCorrV2_NNLO/"
            }
categories = ["CR-mumu"]
ConsiderUnc = True
UseTGraph = True

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

hists = {}
fits = {}
fitunc = {}
Results = {}
varss = ["ptll_fine"]
inputs = {}
for era in inputsdir:
  inputs[era] = []
  for d in os.walk(inputsdir[era]):
    inputs[era] += [os.path.join(d[0],f) for f in d[2] if f.endswith(".root")]

processes = []
for era in inputs:
  for i in [f for f in inputs[era] if f.endswith(f"{categories[0]}.root")]:
    if f"HTauTau_{era.split('_')[0]}_Hlep" not in i:
      proc = i.split("_CR-mumu")[0].split("/")[-1]
    else:
      proc = i.split(f"_HTauTau_{era.split('_')[0]}_Hlep")[0].split("/")[-1]
    if proc not in processes: processes.append(proc)

def DiagonalizeLinearTransformMatrix(matrix):
    assert matrix.shape[0] == matrix.shape[1], "Matrix is not square."
    _, eigen = np.linalg.eig(matrix)
    eigen_inv = np.linalg.inv(eigen)
    diag = np.matmul(np.matmul(eigen_inv, matrix),eigen)
    return diag, eigen, eigen_inv
def CreateVarianceMatrix(std_dev_vector):
    N = std_dev_vector.shape[0]
    var_matrix = np.zeros([N, N])
    for n in range(N):
      var_matrix[n, n] = math.pow(std_dev_vector[n], 2)
    return var_matrix
def MatrixPower(matrix, power):
    diag, eigen, eigen_inv = DiagonalizeLinearTransformMatrix(matrix)
    for n in range(diag.shape[0]):
      diag[n, n] = math.pow(diag[n, n], power)
    return np.matmul(np.matmul(eigen, diag), eigen_inv)
def CovWhitening(cov_matrix):
    N = cov_matrix.shape[0]
    assert N == cov_matrix.shape[1], "Covariance matrix is not square."
    corr_matrix = np.zeros([N, N])
    std_dev_vector = np.zeros([N])
    for n in range(N):
      assert cov_matrix[n, n] > 0, "Variance is not positive."
      std_dev_vector[n] = math.sqrt(cov_matrix[n, n])

    for n in range(N):
      corr_matrix[n][n] = 1
      for k in range(n,N):
        assert cov_matrix[n, k] == cov_matrix[k, n], "Covariance matrix is not symmetric."
        corr_matrix[n, k] = cov_matrix[n, k] / (std_dev_vector[n] * std_dev_vector[k])
        corr_matrix[k, n] = corr_matrix[n, k]

    assert corr_matrix.shape[1] == N, "Correlation matrix is not square."
    assert std_dev_vector.shape[0] == N, "Size of the std dev vector do not correspond to the size of the correlation matrix."

    corr_inv_sqrt = MatrixPower(corr_matrix, -0.5)
    var_matrix = CreateVarianceMatrix(std_dev_vector)
    var_inv_sqrt = MatrixPower(var_matrix, -0.5)
    return np.linalg.inv(np.matmul(corr_inv_sqrt, var_inv_sqrt))

for era in inputsdir:
  acc = era.split("_")[1]
  hists[era] = {}
  for c in categories:
    hists[era][c] = {}
    for f in inputs[era]:
      if not f.endswith(f"{c}.root"): continue
      proc = ""
      for p in processes:
        if f.split("/")[-1].startswith(p):
          proc = p
          break
      rfile = ROOT.TFile(os.path.join(f), "read")
      histnameaddition = proc if ("part" not in f.split("/")[-1]) else f.split(f"_HTauTau_{era.split('_')[0]}_Hlep_CR-mumu")[0].split("/")[-1]+"proot"
      if proc.startswith("Muon"): proc = "Muon" # Combine data eras
      if proc.startswith("DY"): proc = "DY" # Combine DY
      for var in varss:
        h = rfile.Get(f"{var}_{histnameaddition}")
        h.SetDirectory(0)
        if proc not in hists[era][c]:
          hists[era][c][proc]={}
        if var not in hists[era][c][proc]:
          hists[era][c][proc][var] = h
          hists[era][c][proc][var].SetDirectory(0)
          # Basic plotting settings
          hists[era][c][proc][var].SetLineColor(602)
          hists[era][c][proc][var].SetLineWidth(1)
          hists[era][c][proc][var].SetFillColor(0)
          hists[era][c][proc][var].SetMarkerStyle(1)
          hists[era][c][proc][var].SetMarkerSize(1.0)
          hists[era][c][proc][var].GetXaxis().SetTitleSize(0.035)
          hists[era][c][proc][var].GetYaxis().SetTitleSize(0.035)
          hists[era][c][proc][var].GetXaxis().SetTitleOffset(1.0)
          hists[era][c][proc][var].GetYaxis().SetTitleOffset(0.0)
          hists[era][c][proc][var].GetXaxis().SetLabelSize(0.035)
          hists[era][c][proc][var].GetYaxis().SetLabelSize(0.035)
          hists[era][c][proc][var].GetXaxis().SetLabelOffset(0.005)
          hists[era][c][proc][var].GetYaxis().SetLabelOffset(0.005)
        else:
          hists[era][c][proc][var].Add(h)
      rfile.Close()

  # Make DataMinusMC
  for c in hists[era]:
    hists[era][c]["DataMinusMC"] = {}
    for var in varss:
      hists[era][c]["DataMinusMC"][var] = hists[era][c]["Muon"][var].Clone()
    for proc in hists[era][c]:
      if proc in ["DY", "Muon", "DataMinusMC"]: continue
      for var in hists[era][c][proc]:
        hists[era][c]["DataMinusMC"][var].Add(hists[era][c][proc][var], -1)

  # Rebin
  binning = {}
  binning["CR-mumu"] = [i/2 for i in range(0,20,1)] + list(range(10,30,1)) + list(range(30,50,2)) + list(range(50,110,5)) + list(range(110,140,10)) + list(range(140,200,30)) + list(range(200,500,100)) + list(range(500,1000,500))

  binning_array = {}
  for c in binning:
    binning[c] += [hists[era][c]["DataMinusMC"][var].GetXaxis().GetBinUpEdge(hists[era][c]["DataMinusMC"][var].GetNbinsX())]
    binning_array[c] = array("d", binning[c])

  # Rebin and print new bin content and errors
  for c in hists[era]:
    for proc in ["DataMinusMC", "DY"]:
      hists[era][c][proc+"_rebin"] = {}
      for var in varss:
        hists[era][c][proc+"_rebin"][var] = hists[era][c][proc][var].Rebin(len(binning[c])-1, hists[era][c][proc][var].GetName()+"_rebin", binning_array[c])
        if proc=="DY":
          for i in range(hists[era][c][proc+"_rebin"][var].GetNbinsX()):
            err = hists[era][c][proc+"_rebin"][var].GetBinError(i+1)/hists[era][c][proc+"_rebin"][var].GetBinContent(i+1) if hists[era][c][proc+"_rebin"][var].GetBinContent(i+1)>0 else 0

  for c in hists[era]:
    hists[era][c]["Ratio_rebin"] = {}
    for var in varss:
      hists[era][c]["Ratio_rebin"][var] = hists[era][c]["DataMinusMC_rebin"][var].Clone(f"{era}_{c}_{var}_Ratio_rebin")
      hists[era][c]["Ratio_rebin"][var].Divide(hists[era][c]["DY_rebin"][var])

  # Make TGraph to put x value in mean position
  if UseTGraph:
    UseTGraph = "Graph"
    MeanAlongX = {}
    for c in hists[era]:
      MeanAlongX[c] = {}
      for var in varss:
        MeanAlongX[c][var] = []
        for i in range(len(binning[c])-1):
          temp = {}
          summ = 0.0
          lowbin = binning[c][i]
          highbin = binning[c][i+1]
          for proc in ["DataMinusMC", "DY"]:
            temp[proc] = {}
            for j in range(hists[era][c][proc][var].GetNbinsX()):
              cent = hists[era][c][proc][var].GetBinCenter(j+1)
              if cent>lowbin and cent<highbin:
                temp[proc][cent] = hists[era][c][proc][var].GetBinContent(j+1)
                summ += temp[proc][cent]
              elif cent>highbin: break
          val = sum([cent*temp[proc][cent]/summ for cent in temp[proc] for proc in temp])
          MeanAlongX[c][var].append(val)
    for c in hists[era]:
      hists[era][c]["Ratio_rebinGraph"] = {}
      for var in varss:
        xvals = MeanAlongX[c][var]
        assert hists[era][c]["Ratio_rebin"][var].GetNbinsX() == len(MeanAlongX[c][var])
        yvals = [hists[era][c]["Ratio_rebin"][var].GetBinContent(i+1) for i in range(hists[era][c]["Ratio_rebin"][var].GetNbinsX())]
        xerrdown = [MeanAlongX[c][var][i]-binning[c][i] for i in range(len(MeanAlongX[c][var]))]
        xerrup = [binning[c][i+1]-MeanAlongX[c][var][i] for i in range(len(MeanAlongX[c][var]))]
        yerr = [hists[era][c]["Ratio_rebin"][var].GetBinError(i+1) for i in range(hists[era][c]["Ratio_rebin"][var].GetNbinsX())]
        hists[era][c]["Ratio_rebinGraph"][var] = ROOT.TGraphAsymmErrors(len(xvals), array('d', xvals), array('d', yvals), array('d', xerrdown), array('d', xerrup), array('d', yerr), array('d', yerr))
        hists[era][c]["Ratio_rebinGraph"][var].SetTitle(f"{era}_{c}_{var}_Ratio_rebinGraph")

  # Undo normalization factor
  norm = {}
  for c in hists[era]:
    for var in varss:
      norm[c+"_"+var] = hists[era][c]["DY_rebin"][var].Integral() / hists[era][c]["DataMinusMC_rebin"][var].Integral()

  # Make fits
  Results[era] = {}
  fits[era] = {}
  fitunc[era] = {}
  for c in hists[era]:
    Results[era][c] = {}
    fits[era][c] = {}
    fitunc[era][c] = {}
    for var in varss:
      Results[era][c][var] = {}
      fitunc[era][c][var] = {}
      for proc in ["DataMinusMC", "DY", "Ratio", "Error"]:
        if proc!="Error":
          Results[era][c][var][proc] = [hists[era][c][proc+"_rebin"][var].GetBinContent(i+1) for i in range(hists[era][c][proc+"_rebin"][var].GetNbinsX())]
          canv = ROOT.TCanvas(f'c_{era}_{c}_{proc}_{var}', f'c_{era}_{c}_{proc}_{var}', 600, 600)
          canv.cd()
          xname = "p_{T,ll} [GeV]"
          if proc!="Ratio":
            hists[era][c][proc+"_rebin"][var].GetXaxis().SetTitle(xname)
            hists[era][c][proc+"_rebin"][var].GetYaxis().SetTitle("Events")
            title = f"DY {acc} MC" if proc=="DY" else "Data minus non-DY MC"
            hists[era][c][proc+"_rebin"][var].SetTitle(title)
            hists[era][c][proc+"_rebin"][var].Draw("HIST")
          else:
            hists[era][c][proc+"_rebin"+UseTGraph][var].GetXaxis().SetTitle(xname)
            hists[era][c][proc+"_rebin"+UseTGraph][var].GetYaxis().SetTitle("Ratio (Data-MC)/DY")
            hists[era][c][proc+"_rebin"+UseTGraph][var].SetTitle("Ratio")
            if not UseTGraph:
              hists[era][c][proc+"_rebin"+UseTGraph][var].Draw("EP")
            else:
              hists[era][c][proc+"_rebin"+UseTGraph][var].Draw("AP")
          canv.SetTopMargin(0.075)
          canv.SetBottomMargin(0.15)
          canv.SetLeftMargin(0.15)
          canv.SetRightMargin(0.05)
          if not os.path.exists(f"plots_{era}"): os.makedirs(f"plots_{era}")
          if not os.path.exists(f"plots_{era}/hists"): os.makedirs(f"plots_{era}/hists")
          canv.SaveAs(f"plots_{era}/hists/Hist_{era}_{c}_{proc}_{var}.png")
        else:
          Results[era][c][var][proc] = [hists[era][c]["Ratio_rebin"][var].GetBinError(i+1) for i in range(hists[era][c]["Ratio_rebin"][var].GetNbinsX())]
        if proc=="Ratio":
          binning = [hists[era][c][proc+"_rebin"][var].GetBinLowEdge(i+1) for i in range(hists[era][c][proc+"_rebin"][var].GetNbinsX())]
          binning.append(binning[-1] + hists[era][c][proc+"_rebin"][var].GetBinWidth(hists[era][c][proc+"_rebin"][var].GetNbinsX()+1))
          Results[era][c]["HistBinning"] = binning

          if acc=="NLO":
            # Gaus+Gaus+Linear
            fits[era][c][var] = ROOT.TF1(f"fit_{era}_{c}_{var}", "([0]*exp(-0.5*((x-[1])/[2])^2)+[3])*(erf(-999*(x-[8]))+1)/2 + ([0]*exp(-0.5*(([8]-[1])/[2])^2)+[3]+[4]*exp(-0.5*((x-[5])/[6])^2)-[4]*exp(-0.5*(([8]-[5])/[6])^2))*(erf(999*(x-[8]))+1)/2*(erf(-999*(x-[9]))+1)/2 + ([0]*exp(-0.5*(([8]-[1])/[2])^2)+[3]+[4]*exp(-0.5*(([9]-[5])/[6])^2)-[4]*exp(-0.5*(([8]-[5])/[6])^2)+x*[7]+[9]*[7])*(erf(999*(x-[9]))+1)/2", binning[0], binning[-1])
            fits[era][c][var].SetParameter(0,-0.2)
            fits[era][c][var].SetParameter(1,5.0)
            fits[era][c][var].SetParameter(2,3.0)
            fits[era][c][var].SetParameter(3,1.1)
            fits[era][c][var].SetParameter(4,-0.25)
            fits[era][c][var].SetParameter(5,40)
            fits[era][c][var].SetParameter(6,16.0)
            fits[era][c][var].SetParameter(7,-0.0001)
            fits[era][c][var].SetParameter(8,15)
            fits[era][c][var].SetParameter(9,50)
          elif acc=="LO":
            # Gaus+Log+Linear
            fits[era][c][var] = ROOT.TF1(f"fit_{era}_{c}_{var}", "([0]*exp(-0.5*((x-[1])/[2])^2)+[3])*(erf(-999*(x-[5]))+1)/2 + ([4]*log(x)+[0]*exp(-0.5*(([5]-[1])/[2])^2)+[3]-[4]*log([5]))*(erf(999*(x-[5]))+1)/2*(erf(-999*(x-[6]))+1)/2 + ([4]*log([6])+[0]*exp(-0.5*(([5]-[1])/[2])^2)+[3]-[4]*log([5])+[7]*x-[6]*[7])*(erf(999*(x-[6]))+1)/2", binning[0], binning[-1])
            fits[era][c][var].SetParameter(0,-0.32)
            fits[era][c][var].SetParameter(1,1.2)
            fits[era][c][var].SetParameter(2,3.8)
            fits[era][c][var].SetParameter(3,1.0)
            fits[era][c][var].SetParameter(4,0.15)
            fits[era][c][var].SetParameter(5,13)
            fits[era][c][var].SetParameter(6,120)
            fits[era][c][var].SetParameter(7,-0.00059)
          elif acc=="NNLO":
            # GausLin+Linear+Log
            fits[era][c][var] = ROOT.TF1(f"fit_{era}_{c}_{var}", "([0]*exp(-0.5*((x-[1])/[2])^2)+[3]*x+[4])*(erf(-999*(x-[7]))+1)/2 + ([0]*exp(-0.5*(([7]-[1])/[2])^2)+[3]*[7]+[4]+[5]*x-[5]*[7])*(erf(999*(x-[7]))+1)/2*(erf(-999*(x-[8]))+1)/2 + ([0]*exp(-0.5*(([7]-[1])/[2])^2)+[3]*[7]+[4]+[5]*[8]-[5]*[7]+[6]*log(x)-[6]*log([8]))*(erf(999*(x-[8]))+1)/2", binning[0], binning[-1])
            fits[era][c][var].SetParameter(0,-0.2)
            fits[era][c][var].SetParameter(1,5.0)
            fits[era][c][var].SetParameter(2,3.0)
            fits[era][c][var].SetParameter(3,-0.01)
            fits[era][c][var].SetParameter(4,1.1)
            fits[era][c][var].SetParameter(5,0.0)
            fits[era][c][var].SetParameter(6,0.15)
            fits[era][c][var].SetParameter(7,10)
            fits[era][c][var].SetParameter(8,60)
          #ROOT.Math.MinimizerOptions().PrintDefault()
          #ROOT.Math.MinimizerOptions().SetDefaultStrategy(2)
          toll = 0.01
          attempt = 0
          while True:
            #ROOT.Math.MinimizerOptions().SetDefaultMinimizer("Minuit2")
            ROOT.Math.MinimizerOptions().SetDefaultTolerance(toll)
            ROOT.Math.MinimizerOptions().SetDefaultMaxFunctionCalls(1000000)
            ROOT.Math.MinimizerOptions().SetDefaultMaxIterations(1000000)
            ROOT.Math.MinimizerOptions().PrintDefault()
            res = hists[era][c]["Ratio_rebin"+UseTGraph][var].Fit(f"fit_{era}_{c}_{var}", "S0")
            if res.Status()==0: break
            attempt += 1
            toll *= math.sqrt(math.sqrt(10))
            if attempt==20: raise RuntimeError
          for i in range(fits[era][c][var].GetNpar()):
          if ConsiderUnc:
            cor = res.GetCorrelationMatrix()
            cov = res.GetCovarianceMatrix()
            npcov = np.matrix([[cov[x][y] for x in range(fits[era][c][var].GetNpar())] for y in range(fits[era][c][var].GetNpar())])
            white = CovWhitening(npcov)
            for x in range(fits[era][c][var].GetNpar()):
              for y in range(fits[era][c][var].GetNpar()):
                cov[x][y] = white[x, y]
          for i in range(-1,fits[era][c][var].GetNpar()):
            addition = ""
            titleadd = ""
            if i>=0 and not ConsiderUnc: break
            if i>=0:
              addition = f"_unc{i}"
              titleadd = f", Uncertainty {i+1}"
              for unc in ["up","down"]:
                fitunc[era][c][var][f"{i}_{unc}"] = fits[era][c][var].Clone(fits[era][c][var].GetName()+f"_Unc_{i}_{unc}")
                for j in range(fits[era][c][var].GetNpar()):
                  if unc=="up":
                    fitunc[era][c][var][f"{i}_{unc}"].FixParameter(j, fits[era][c][var].GetParameter(j)+cov[j][i])
                  else:
                    fitunc[era][c][var][f"{i}_{unc}"].FixParameter(j, fits[era][c][var].GetParameter(j)-cov[j][i])
                fitunc[era][c][var][f"{i}_{unc}"].SetLineColor(3)
                fitunc[era][c][var][f"{i}_{unc}"].SetLineWidth(1)
            canv = ROOT.TCanvas(f'fit_{era}_{c}_{proc}_{var}{addition}', f'fit{era}_{c}_{proc}_{var}{addition}', 600, 600)
            canv.cd()
            hists[era][c]["Ratio_rebin"+UseTGraph][var].SetTitle(f"Fit over weight for {era}"+titleadd)
            hists[era][c]["Ratio_rebin"+UseTGraph][var].GetXaxis().SetRangeUser(2, binning[-1])
            if not UseTGraph:
              hists[era][c]["Ratio_rebin"+UseTGraph][var].Draw("EP")
            else:
              hists[era][c]["Ratio_rebin"+UseTGraph][var].Draw("AP")
            canv.SetTopMargin(0.075)
            canv.SetBottomMargin(0.15)
            canv.SetLeftMargin(0.15)
            canv.SetRightMargin(0.05)
            fits[era][c][var].Draw("SAME")
            if i>=0:
              for unc in ["up","down"]:
                fitunc[era][c][var][f"{i}_{unc}"].Draw("SAME")
            if not os.path.exists(f"plots_{era}/fits"): os.makedirs(f"plots_{era}/fits")
            canv.SaveAs(f"plots_{era}/fits/Fit_{era}_{c}_{proc}_{var}{addition}.png")
            canv.SetLogx()
            if not os.path.exists(f"plots_{era}/fits_log"): os.makedirs(f"plots_{era}/fits_log")
            canv.SaveAs(f"plots_{era}/fits_log/Fit_{era}_{c}_{proc}_{var}{addition}_log.png")

with open('pTllCorrResults.json', 'w') as file:
  json.dump(Results, file)

file = ROOT.TFile.Open("pTllCorrResults.root", "RECREATE")
for era in fits:
  for c in fits[era]:
    for var in fits[era][c]:
      file.WriteObject(fits[era][c][var], "Fit_"+era+"_"+c+"_"+var)
      for unc in fitunc[era][c][var]:
        file.WriteObject(fitunc[era][c][var][unc], "FitUncertainty_"+unc+"_"+era+"_"+c+"_"+var)
file.Close()
exit()

