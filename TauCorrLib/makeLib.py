import os
import json
import ROOT
import math
import re
from array import array
import copy

correctionlibTemplate = {}
correctionlibTemplate["schema_version"] = 2
correctionlibTemplate["corrections"] = [] # list of dicts
correctionlibs = {"Full" : {}, "DY" : {}, "REC" : {}}
for typ in correctionlibs:
  correctionlibs[typ]["AllEras"] = copy.deepcopy(correctionlibTemplate)

def SaveCorrLib(corrlib, outname):
  with open(outname, "w") as outfile:
    json.dump(corrlib, outfile, indent=1)
  os.system(f"gzip {outname}")
  with open(outname, "w") as outfile:
    json.dump(corrlib, outfile, indent=1)

################################################## PTLL INPUTS

ijson = "../DYpTllCorr/pTllCorrResults.json"
with open(ijson) as json_file:
  ptlldata = json.load(json_file)
iroot = "../DYpTllCorr/pTllCorrResults.root"
rfile = ROOT.TFile(iroot)

keynames = [e.GetName() for e in rfile.GetListOfKeys()]
fits = {}
def GetFormula(tf1):
  res = str(tf1.GetExpFormula())
  n = tf1.GetNpar()
  for i in range(n):
    res = res.replace(f"[p{i}]",str(tf1.GetParameter(i)))
  res = res.replace("+-","-").replace("--","+")
  return res
for era in ptlldata:
  fits[era] = {}
  for fs in ptlldata[era]:
    fits[era][fs] = {}
    for var in ptlldata[era][fs]:
      if var=="HistBinning": continue
      fits[era][fs][var] = {}
      fits[era][fs][var]["nom"] = GetFormula(rfile.Get(f"Fit_{era}_{fs}_{var}"))
      uncs = [k for k in keynames if k.startswith("FitUncertainty") and k.endswith(f"_{era}_{fs}_{var}")]
      for unc in uncs:
        pattern = re.match(f"FitUncertainty_([0-9]+)_(up|down)_{era}_{fs}_{var}", unc)
        fits[era][fs][var][pattern.group(2)+str(int(pattern.group(1))+1)] = GetFormula(rfile.Get(unc))
      fits[era][fs][var]["N_unc"]=int(len(uncs)/2)

################################################## PTLL CORRLIB

alleras = list(dict.fromkeys([f'{era.split("_")[0]}' for era in ptlldata]))
allorders = sorted(list(dict.fromkeys([f'{era.split("_")[1]}' for era in ptlldata])))
for typ in ["Full", "DY"]:
  for era in alleras:
    if era not in correctionlibs[typ]:
      correctionlibs[typ][era] = copy.deepcopy(correctionlibTemplate)
PTLL = {}
PTLL["name"] = "DY_pTll_reweighting"
PTLL["description"] = "Weights that are to be applied on DY MC samples. The weights were derived from a mu-mu control region."
PTLL["version"] = 1
PTLL["inputs"] = []
PTLL["inputs"].append({"name": "era", "type": "string", "description": f"Era: {', '.join(alleras)}"})
PTLL["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
PTLL["inputs"].append({"name": "ptll", "type": "real", "description": "Gen-level pTll: Obtained from gen-level Electrons/Muons with status==1, or gen-level Taus with status==2, both also requiring status flag 'fromHardProcess'"})
PTLL["inputs"].append({"name": "syst", "type": "string", "description": "Systematic variations: 'nom', 'up1', 'down1', 'up2', 'down2', ..."})
PTLL["output"] = {"name": "sf", "type": "real", "description": "Resulting weight"}
PTLL["data"] = {}
PTLL["data"]["nodetype"] = "category"
PTLL["data"]["input"] = "era"
PTLL["data"]["content"] = []
plain = copy.deepcopy(PTLL)

for actera in alleras:
 ptllera = {}
 ptllera["key"] = actera
 ptllera["value"] = {}
 ptllera["value"]["nodetype"] = "category"
 ptllera["value"]["input"] = "order"
 ptllera["value"]["content"] = []
 byeraPTLL = copy.deepcopy(plain)
 del byeraPTLL["inputs"][0]
 byeraPTLL["data"]["input"] = "order"
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in ptlldata: continue
  orderbin = {}
  orderbin["key"] = order
  orderbin["value"] = {}

  # Weight from fit
  orderbin["value"]["nodetype"] = "category"
  orderbin["value"]["input"] = "syst"
  orderbin["value"]["content"] = []
  for key in ["nom"]+[f"{updn}{i+1}" for updn in ["up", "down"] for i in range(fits[era]["CR-mumu"]["ptll_fine"]["N_unc"])]:
    sysbin = {}
    sysbin["key"] = key
    sysbin["value"] = {}
    sysbin["value"]["nodetype"] = "formula"
    sysbin["value"]["expression"] = fits[era]["CR-mumu"]["ptll_fine"][key]
    sysbin["value"]["parser"] = "TFormula"
    sysbin["value"]["variables"] = ["ptll"]
    orderbin["value"]["content"].append(sysbin)
  ptllera["value"]["content"].append(orderbin)
  byeraPTLL["data"]["content"].append(orderbin)
 correctionlibs["Full"][actera]["corrections"].append(byeraPTLL)
 correctionlibs["DY"][actera]["corrections"].append(byeraPTLL)
 PTLL["data"]["content"].append(ptllera)

correctionlibs["Full"]["AllEras"]["corrections"].append(PTLL)
correctionlibs["DY"]["AllEras"]["corrections"].append(PTLL)

################################################## PTLL CORRLIB N_UNC

PTLL = {}
PTLL["name"] = "DY_pTll_reweighting_N_uncertainty"
PTLL["description"] = "Just gives the number of uncertainties (varies per DY samples type: 'LO'=madgraph, 'NLO'=amcatnlo, 'NNLO'=powheg). Output must be converted to integer."
PTLL["version"] = 1
PTLL["inputs"] = []
PTLL["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
PTLL["output"] = {"name": "nunc", "type": "real", "description": "Number of uncertainties"}
PTLL["data"] = {}
PTLL["data"]["nodetype"] = "category"
PTLL["data"]["input"] = "order"
PTLL["data"]["content"] = []

for actera in alleras:
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in ptlldata: continue
  if order in [alpha["key"] for alpha in PTLL["data"]["content"]]: continue # Same for all eras; when one order done, just skip
  ptllera = {}
  ptllera["key"] = order
  ptllera["value"] = float(fits[era]["CR-mumu"]["ptll_fine"]["N_unc"])
  PTLL["data"]["content"].append(ptllera)
 correctionlibs["Full"][actera]["corrections"].append(PTLL)
 correctionlibs["DY"][actera]["corrections"].append(PTLL)

correctionlibs["Full"]["AllEras"]["corrections"].append(PTLL)
correctionlibs["DY"]["AllEras"]["corrections"].append(PTLL)


################################################## RECOIL CORR INPUTS

ijson = "../RecoilCorr/RecoilCorrFitResults.json"
with open(ijson) as json_file:
  recoildata = json.load(json_file)
iroot = "../RecoilCorr/RecoilCorrFitResults.root"
rfile = ROOT.TFile(iroot)

hists = {}
for fs in recoildata:
  hists[fs] = {}
  for era in recoildata[fs]:
    hists[fs][era] = {}
    for category in recoildata[fs][era]:
      if category=="HistBinning": continue
      hists[fs][era][category] = {}
      for kind in ["upara", "uperp"]:
        hists[fs][era][category][kind] = {}
        hists[fs][era][category][kind]["Data"] = rfile.Get(f"Hist_{fs}_{era}_{category}_{kind}_DataMinusMC")
        hists[fs][era][category][kind]["DY"] = rfile.Get(f"Hist_{fs}_{era}_{category}_{kind}_DY")

fits = {}
unc = {}
fs = "mumu"
for era in recoildata[fs]:
  fits[era] = {}
  unc[era] = {}
  unc[era]["Response"] = {}
  for category in recoildata[fs][era]:
    fits[era][category] = {}
    if "Hist" in category: continue
    if "ZpT" not in category: continue
    for kind in ["upara", "uperp"]:
      fits[era][category][kind] = {}
      for data in ["Data", "DY"]:
        fits[era][category][kind][data] = {}
        readdata = "DataMinusMC" if data=="Data" else "DY"
        fits[era][category][kind][data]["fit"] = rfile.Get(f"Fit_{fs}_{era}_{category}_{kind}_{readdata}")
        if kind == "uperp":
          fits[era][category][kind][data]["cdf"] = ROOT.TF1(f"CDF_{era}_{category}_{kind}_{readdata}","[0]/2*(1+erf(x/[1]))+[2]/2*(1+erf(x/[3]))", fits[era][category][kind][data]["fit"].GetXmin(), fits[era][category][kind][data]["fit"].GetXmax());
          sigma = fits[era][category][kind][data]["fit"].GetParameter(1)
          fits[era][category][kind][data]["cdf"].FixParameter(0, fits[era][category][kind][data]["fit"].GetParameter(0)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(1, sigma*math.sqrt(2))
          sigma = fits[era][category][kind][data]["fit"].GetParameter(3)
          fits[era][category][kind][data]["cdf"].FixParameter(2, fits[era][category][kind][data]["fit"].GetParameter(2)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(3, sigma*math.sqrt(2))
        elif kind=="upara":
          fits[era][category][kind][data]["cdf"] = ROOT.TF1(f"CDF_{era}_{category}_{kind}_{readdata}","[0]/2*(1+erf((x-[1])/[2]))*(x<[1])+([3]/2*(1+erf((x-[1])/[4]))+([0]-[3])/2)*(x>=[1])+[5]/2*(1+erf((x-[6])/[7]))*(x<[6])+([8]/2*(1+erf((x-[6])/[9]))+([5]-[8])/2)*(x>=[6])", fits[era][category][kind][data]["fit"].GetXmin(), fits[era][category][kind][data]["fit"].GetXmax());
          sigma = fits[era][category][kind][data]["fit"].GetParameter(2)
          fits[era][category][kind][data]["cdf"].FixParameter(0, fits[era][category][kind][data]["fit"].GetParameter(0)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(1, fits[era][category][kind][data]["fit"].GetParameter(1))
          fits[era][category][kind][data]["cdf"].FixParameter(2, sigma*math.sqrt(2))
          sigma = fits[era][category][kind][data]["fit"].GetParameter(3)
          fits[era][category][kind][data]["cdf"].FixParameter(3, fits[era][category][kind][data]["fit"].GetParameter(0)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(4, sigma*math.sqrt(2))
          sigma = fits[era][category][kind][data]["fit"].GetParameter(6)
          fits[era][category][kind][data]["cdf"].FixParameter(5, fits[era][category][kind][data]["fit"].GetParameter(4)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(6, fits[era][category][kind][data]["fit"].GetParameter(5))
          fits[era][category][kind][data]["cdf"].FixParameter(7, sigma*math.sqrt(2))
          sigma = fits[era][category][kind][data]["fit"].GetParameter(7)
          fits[era][category][kind][data]["cdf"].FixParameter(8, fits[era][category][kind][data]["fit"].GetParameter(4)*math.sqrt(2*math.pi)*sigma)
          fits[era][category][kind][data]["cdf"].FixParameter(9, sigma*math.sqrt(2))
    if "response" in recoildata[fs][era][category]: unc[era]["Response"][category] = recoildata[fs][era][category]["response"]["DY"]["ResponseMean"]

for era in recoildata["ee"]:
  unc[era]["Uncertainty"] = {}
  for category in recoildata["ee"][era]: # Should be just 0j, 1j, 2j
    unc[era]["Uncertainty"][category] = {}
    unc[era]["Uncertainty"][category]["response"] = recoildata["ee"][era][category]["response"]["DY"]["ResponseUnc"]
    unc[era]["Uncertainty"][category]["resolution"] = recoildata["ee"][era][category]["resolution"]["DY"]["ResolutionUnc"]
    unc[era]["Response"][category+"_0-5"] = f'({unc[era]["Response"][category+"_ZpT10"]})'
    unc[era]["Response"][category+"_5-15"] = f'({(unc[era]["Response"][category+"_ZpT20"]-unc[era]["Response"][category+"_ZpT10"])/10.0}*y+{unc[era]["Response"][category+"_ZpT10"]-(unc[era]["Response"][category+"_ZpT20"]-unc[era]["Response"][category+"_ZpT10"])/10.0*5.0})'
    unc[era]["Response"][category+"_15-25"] = f'({(unc[era]["Response"][category+"_ZpT30"]-unc[era]["Response"][category+"_ZpT20"])/10.0}*y+{unc[era]["Response"][category+"_ZpT20"]-(unc[era]["Response"][category+"_ZpT30"]-unc[era]["Response"][category+"_ZpT20"])/10.0*15.0})'
    unc[era]["Response"][category+"_25-40"] = f'({(unc[era]["Response"][category+"_ZpT50"]-unc[era]["Response"][category+"_ZpT30"])/15.0}*y+{unc[era]["Response"][category+"_ZpT30"]-(unc[era]["Response"][category+"_ZpT50"]-unc[era]["Response"][category+"_ZpT30"])/15.0*25.0})'
    unc[era]["Response"][category+"_40-75"] = f'({(unc[era]["Response"][category+"_ZpTInf"]-unc[era]["Response"][category+"_ZpT50"])/35.0}*y+{unc[era]["Response"][category+"_ZpT50"]-(unc[era]["Response"][category+"_ZpTInf"]-unc[era]["Response"][category+"_ZpT50"])/35.0*40.0})'
    unc[era]["Response"][category+"_75-Inf"] = f'({unc[era]["Response"][category+"_ZpTInf"]})'
    for c in unc[era]["Response"]:
      if "ZpT" not in c:
        unc[era]["Response"][c] = unc[era]["Response"][c].replace("+-","-")


manual = {}
requiredbinning = {}
cdfbinedges = {}
finalweight = {}
for era in hists[fs]:
  manual[era] = {}
  requiredbinning[era] = {}
  cdfbinedges[era] = {}
  finalweight[era] = {}
  for category in hists[fs][era]:
    manual[era][category] = {}
    requiredbinning[era][category] = {}
    cdfbinedges[era][category] = {}
    finalweight[era][category] = {}
    for kind in hists[fs][era][category]:
      manual[era][category][kind] = {}
      cdfbinedges[era][category][kind] = {}
      finalweight[era][category][kind] = {}
      for data in hists[fs][era][category][kind]:
        manual[era][category][kind][data] = {}
        cdfbinedges[era][category][kind][data] = []
        nbins = hists[fs][era][category][kind][data].GetNbinsX()
        binedges = []
        for i in range(nbins):
          binedges.append(hists[fs][era][category][kind][data].GetBinLowEdge(i+1))
        binedges.append(binedges[-1]+hists[fs][era][category][kind][data].GetBinWidth(nbins))
        integral = hists[fs][era][category][kind][data].Integral()
        underflow = hists[fs][era][category][kind][data].GetBinContent(0)
        if underflow < 0: underflow = 0.0
        integral += underflow
        overflow = hists[fs][era][category][kind][data].GetBinContent(nbins+1)
        if overflow < 0: overflow = 0.0
        integral += overflow
        manual[era][category][kind][data]["cdf"] = {}
        manual[era][category][kind][data]["inv"] = {}
        prev = underflow / integral
        for i in range(nbins):
          xmin = binedges[i]
          xmax = binedges[i+1]
          ymin = prev
          ymax = prev+hists[fs][era][category][kind][data].GetBinContent(i+1)/integral
          assert ymax-1 < 1e-14 , f"ymax = {ymax}"
          if ymax>1: ymax=1.0
          manual[era][category][kind][data]["cdf"][f"{xmin}_{xmax}"] = f"{(ymax-ymin)/(xmax-xmin)}*x+{ymin-(ymax-ymin)/(xmax-xmin)*xmin}".replace("+-","-")
          if (ymax-ymin)!=0.0:
            manual[era][category][kind][data]["inv"][f"{ymin}_{ymax}"] = f"{(xmax-xmin)/(ymax-ymin)}*x+{xmin-(xmax-xmin)/(ymax-ymin)*ymin}".replace("+-","-")
          cdfbinedges[era][category][kind][data].append(f"{ymin}_{ymax}")
          prev = ymax
        assert prev<=1.0

      listofubins = sorted([a for a in manual[era][category][kind]["DY"]["cdf"]], key=lambda x: float(x.split("_")[0]))
      listofDYcdfbins = cdfbinedges[era][category][kind]["DY"]
      listofDATAcdfbins = cdfbinedges[era][category][kind]["Data"]
      requiredbinning[era][category][kind] = [float("-inf")]
      prevbin = float(listofubins[0].split('_')[0])
      requiredbinning[era][category][kind].append(prevbin)
      finalweight[era][category][kind][f"-inf_{prevbin}"] = f'{hists[fs][era][category][kind]["Data"].GetMean()}+(x-{hists[fs][era][category][kind]["DY"].GetMean()})*{hists[fs][era][category][kind]["Data"].GetRMS()}/{hists[fs][era][category][kind]["DY"].GetRMS()}'.replace("--","+")
      j = -1
      def IncreaseJ(j):
        j+=1
        l = float(listofDATAcdfbins[j].split("_")[0])
        h = float(listofDATAcdfbins[j].split("_")[1])
        if l==h: return IncreaseJ(j)
        return j,l,h
      j,lowDATAcdf,highDATAcdf = IncreaseJ(j)
      getout = False
      for i,cdfbin in enumerate(listofubins):
        lowbin = float(cdfbin.split("_")[0])
        highbin = float(cdfbin.split("_")[1])
        lowDYcdf = float(listofDYcdfbins[i].split("_")[0])
        highDYcdf = float(listofDYcdfbins[i].split("_")[1])
        if highDYcdf < lowDATAcdf:
          finalweight[era][category][kind][f"-inf_{highbin}"] = finalweight[era][category][kind][f"-inf_{prevbin}"]
          del finalweight[era][category][kind][f"-inf_{prevbin}"]
          del requiredbinning[era][category][kind][-1]
          requiredbinning[era][category][kind].append(highbin)
          prevbin = highbin
          continue
        while lowDYcdf > highDATAcdf:
          assert i==0
          j,lowDATAcdf,highDATAcdf = IncreaseJ(j)
        while highDATAcdf < highDYcdf:
          assert lowDYcdf <= highDATAcdf, f"{lowDYcdf} <= {highDATAcdf}"
          assert highDYcdf >= highDATAcdf, f"{highDYcdf} >= {highDATAcdf}"
          addthis = float(eval(manual[era][category][kind]["DY"]["inv"][listofDYcdfbins[i]].replace("x",str(highDATAcdf))))
          if addthis not in requiredbinning[era][category][kind]:
            requiredbinning[era][category][kind].append(addthis)
            finalweight[era][category][kind][f"{prevbin}_{addthis}"] = manual[era][category][kind]["Data"]["inv"][listofDATAcdfbins[j]].replace("x","("+manual[era][category][kind]["DY"]["cdf"][cdfbin]+")")
            prevbin = addthis
          j+=1
          if j==len(listofDATAcdfbins):
            getout = True
            break
          lowDATAcdf = float(listofDATAcdfbins[j].split("_")[0])
          highDATAcdf = float(listofDATAcdfbins[j].split("_")[1])
        if getout: break
        if highbin not in requiredbinning[era][category][kind]:
          requiredbinning[era][category][kind].append(highbin)
          finalweight[era][category][kind][f"{prevbin}_{highbin}"] = manual[era][category][kind]["Data"]["inv"][listofDATAcdfbins[j]].replace("x","("+manual[era][category][kind]["DY"]["cdf"][cdfbin]+")")
          prevbin = highbin
      requiredbinning[era][category][kind].append(float("inf"))
      finalweight[era][category][kind][f"{prevbin}_inf"] = f'{hists[fs][era][category][kind]["Data"].GetMean()}+(x-{hists[fs][era][category][kind]["DY"].GetMean()})*{hists[fs][era][category][kind]["Data"].GetRMS()}/{hists[fs][era][category][kind]["DY"].GetRMS()}'.replace("--","+")
      assert requiredbinning[era][category][kind] == sorted(requiredbinning[era][category][kind])

################################################## RECOIL CORR RESCALING

alleras = list(dict.fromkeys([f'{era.split("_")[0]}' for era in recoildata['mumu']]))
allorders = sorted(list(dict.fromkeys([f'{era.split("_")[1]}' for era in recoildata['mumu']])))
for typ in ["Full", "REC"]:
  for era in alleras:
    if era not in correctionlibs[typ]:
      correctionlibs[typ][era] = copy.deepcopy(correctionlibTemplate)
REC = {}
REC["name"] = "Recoil_correction_Rescaling"
REC["description"] = "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_Rescaling' gives Upara/Uperp using Mean/RMS of the distributions."
REC["version"] = 1
REC["inputs"] = []
REC["inputs"].append({"name": "era", "type": "string", "description": f"Era: {', '.join(alleras)}"})
REC["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
REC["inputs"].append({"name": "njet", "type": "real", "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"})
REC["inputs"].append({"name": "ptll", "type": "real", "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"})
REC["inputs"].append({"name": "var", "type": "string", "description": "The variable name you are giving the input for: 'Upara', 'Uperp' (string). The output will be for the same kind of variable."})
REC["inputs"].append({"name": "val", "type": "real", "description": "Input value of either Upara or Uperp"})
REC["output"] = {"name": "outval", "type": "real", "description": "Output value of either Upara or Uperp"}
REC["data"] = {}
REC["data"]["nodetype"] = "category"
REC["data"]["input"] = "era"
REC["data"]["content"] = []
plain = copy.deepcopy(REC)

for actera in alleras:
 recera = {}
 recera["key"] = actera
 recera["value"] = {}
 recera["value"]["nodetype"] = "category"
 recera["value"]["input"] = "order"
 recera["value"]["content"] = []
 byeraREC = copy.deepcopy(plain)
 del byeraREC["inputs"][0]
 byeraREC["data"]["input"] = "order"
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in recoildata["mumu"]: continue
  orderbin = {}
  orderbin["key"] = order
  orderbin["value"] = {}
  orderbin["value"]["nodetype"] = "binning"
  orderbin["value"]["flow"] = "error"
  orderbin["value"]["input"] = "njet"
  orderbin["value"]["edges"] = [-0.5, 0.5, 1.5, float("inf")]
  orderbin["value"]["content"] = []
  for njet in ["0j", "1j", "2j"]:
    jetbin = {}
    jetbin["nodetype"] = "binning"
    jetbin["flow"] = "error"
    jetbin["input"] = "ptll"
    jetbin["edges"] = [0.0, 10.0, 20.0, 30.0, 50.0, float("inf")]
    jetbin["content"] = []
    for ptll in ["ZpT10", "ZpT20", "ZpT30", "ZpT50", "ZpTInf"]:
      ptllbin = {}
      ptllbin["nodetype"] = "category"
      ptllbin["input"] = "var"
      ptllbin["content"] = []
      for var in ["Upara", "Uperp"]:
        varbin = {}
        varbin["key"] = var
        varbin["value"] = {}
        varbin["value"]["nodetype"] = "formula"
        varbin["value"]["expression"] = f'{hists["mumu"][era][njet+"_"+ptll][var.lower()]["Data"].GetMean()}+(x-{hists["mumu"][era][njet+"_"+ptll][var.lower()]["DY"].GetMean()})*{hists["mumu"][era][njet+"_"+ptll][var.lower()]["Data"].GetRMS()}/{hists["mumu"][era][njet+"_"+ptll][var.lower()]["DY"].GetRMS()}'
        varbin["value"]["expression"] = varbin["value"]["expression"].replace("+-","-").replace("--","+")
        varbin["value"]["parser"] = "TFormula"
        varbin["value"]["variables"] = ["val"]
        ptllbin["content"].append(varbin)
      jetbin["content"].append(ptllbin)
    orderbin["value"]["content"].append(jetbin)
  recera["value"]["content"].append(orderbin)
  byeraREC["data"]["content"].append(orderbin)
 correctionlibs["Full"][actera]["corrections"].append(byeraREC)
 correctionlibs["REC"][actera]["corrections"].append(byeraREC)
 REC["data"]["content"].append(recera)

correctionlibs["Full"]["AllEras"]["corrections"].append(REC)
correctionlibs["REC"]["AllEras"]["corrections"].append(REC)

################################################## RECOIL CORR QUANTILE FIT

REC = {}
REC["name"] = "Recoil_correction_QuantileMapFit"
REC["description"] = "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_QuantileMapFit' returns CDF value from the DY/Data fit (NOT new Upara/Uperp values!) as function of input value. The idea is to get the CDF value from the DY distribution, then find the new Upara/Uperp value that gives the same CDF value from the Data distribution"
REC["version"] = 1
REC["inputs"] = []
REC["inputs"].append({"name": "era", "type": "string", "description": f"Era: {', '.join(alleras)}"})
REC["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
REC["inputs"].append({"name": "njet", "type": "real", "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"})
REC["inputs"].append({"name": "ptll", "type": "real", "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"})
REC["inputs"].append({"name": "var", "type": "string", "description": "The variable name you are giving the input for: 'Upara', 'Uperp' (string)"})
REC["inputs"].append({"name": "kind", "type": "string", "description": "Get CDF value from either 'Data' or 'DY' distribution"})
REC["inputs"].append({"name": "val", "type": "real", "description": "Input value of either Upara or Uperp"})
REC["output"] = {"name": "outval", "type": "real", "description": "Output value from CDF"}
REC["data"] = {}
REC["data"]["nodetype"] = "category"
REC["data"]["input"] = "era"
REC["data"]["content"] = []
plain = copy.deepcopy(REC)

for actera in alleras:
 recera = {}
 recera["key"] = actera
 recera["value"] = {}
 recera["value"]["nodetype"] = "category"
 recera["value"]["input"] = "order"
 recera["value"]["content"] = []
 byeraREC = copy.deepcopy(plain)
 del byeraREC["inputs"][0]
 byeraREC["data"]["input"] = "order"
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in recoildata["mumu"]: continue
  orderbin = {}
  orderbin["key"] = order
  orderbin["value"] = {}
  orderbin["value"]["nodetype"] = "binning"
  orderbin["value"]["flow"] = "error"
  orderbin["value"]["input"] = "njet"
  orderbin["value"]["edges"] = [-0.5, 0.5, 1.5, float("inf")]
  orderbin["value"]["content"] = []
  for njet in ["0j", "1j", "2j"]:
    jetbin = {}
    jetbin["nodetype"] = "binning"
    jetbin["flow"] = "error"
    jetbin["input"] = "ptll"
    jetbin["edges"] = [0.0, 10.0, 20.0, 30.0, 50.0, float("inf")]
    jetbin["content"] = []
    for ptll in ["ZpT10", "ZpT20", "ZpT30", "ZpT50", "ZpTInf"]:
      ptllbin = {}
      ptllbin["nodetype"] = "category"
      ptllbin["input"] = "var"
      ptllbin["content"] = []
      for var in ["Upara", "Uperp"]:
        varbin = {}
        varbin["key"] = var
        varbin["value"] = {}
        varbin["value"]["nodetype"] = "category"
        varbin["value"]["input"] = "kind"
        varbin["value"]["content"] = []
        for kind in ["DY", "Data"]:
          kindbin = {}
          kindbin["key"] = kind
          kindbin["value"] = {}
          kindbin["value"]["nodetype"] = "formula"
          kindbin["value"]["expression"] = GetFormula(fits[era][njet+"_"+ptll][var.lower()][kind]["cdf"])
          kindbin["value"]["parser"] = "TFormula"
          kindbin["value"]["variables"] = ["val"]
          varbin["value"]["content"].append(kindbin)
        ptllbin["content"].append(varbin)
      jetbin["content"].append(ptllbin)
    orderbin["value"]["content"].append(jetbin)
  recera["value"]["content"].append(orderbin)
  byeraREC["data"]["content"].append(orderbin)
 correctionlibs["Full"][actera]["corrections"].append(byeraREC)
 correctionlibs["REC"][actera]["corrections"].append(byeraREC)
 REC["data"]["content"].append(recera)

correctionlibs["Full"]["AllEras"]["corrections"].append(REC)
correctionlibs["REC"]["AllEras"]["corrections"].append(REC)

################################################## RECOIL CORR QUANTILE HIST

REC = {}
REC["name"] = "Recoil_correction_QuantileMapHist"
REC["description"] = "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_QuantileMapHist' returns Upara/Uperp using the CDF obtained directly from histogram instead."
REC["version"] = 1
REC["inputs"] = []
REC["inputs"].append({"name": "era", "type": "string", "description": f"Era: {', '.join(alleras)}"})
REC["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
REC["inputs"].append({"name": "njet", "type": "real", "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"})
REC["inputs"].append({"name": "ptll", "type": "real", "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"})
REC["inputs"].append({"name": "var", "type": "string", "description": "The variable name you are giving the input for: 'Upara', 'Uperp' (string). The output will be for the same kind of variable."})
REC["inputs"].append({"name": "val", "type": "real", "description": "Input value of either Upara or Uperp"})
REC["output"] = {"name": "outval", "type": "real", "description": "Output value of either Upara or Uperp"}
REC["data"] = {}
REC["data"]["nodetype"] = "category"
REC["data"]["input"] = "era"
REC["data"]["content"] = []
plain = copy.deepcopy(REC)

for actera in alleras:
 recera = {}
 recera["key"] = actera
 recera["value"] = {}
 recera["value"]["nodetype"] = "category"
 recera["value"]["input"] = "order"
 recera["value"]["content"] = []
 byeraREC = copy.deepcopy(plain)
 del byeraREC["inputs"][0]
 byeraREC["data"]["input"] = "order"
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in recoildata["mumu"]: continue
  orderbin = {}
  orderbin["key"] = order
  orderbin["value"] = {}
  orderbin["value"]["nodetype"] = "binning"
  orderbin["value"]["flow"] = "error"
  orderbin["value"]["input"] = "njet"
  orderbin["value"]["edges"] = [-0.5, 0.5, 1.5, float("inf")]
  orderbin["value"]["content"] = []
  for njet in ["0j", "1j", "2j"]:
    jetbin = {}
    jetbin["nodetype"] = "binning"
    jetbin["flow"] = "error"
    jetbin["input"] = "ptll"
    jetbin["edges"] = [0.0, 10.0, 20.0, 30.0, 50.0, float("inf")]
    jetbin["content"] = []
    for ptll in ["ZpT10", "ZpT20", "ZpT30", "ZpT50", "ZpTInf"]:
      ptllbin = {}
      ptllbin["nodetype"] = "category"
      ptllbin["input"] = "var"
      ptllbin["content"] = []
      for var in ["Upara", "Uperp"]:
        varbin = {}
        varbin["key"] = var
        varbin["value"] = {}
        varbin["value"]["nodetype"] = "binning"
        varbin["value"]["flow"] = "error"
        varbin["value"]["input"] = "val"
        varbin["value"]["edges"] = requiredbinning[era][njet+"_"+ptll][var.lower()]
        varbin["value"]["content"] = []
        assert len(requiredbinning[era][njet+"_"+ptll][var.lower()]) == len(finalweight[era][njet+"_"+ptll][var.lower()])+1
        for i in range(len(requiredbinning[era][njet+"_"+ptll][var.lower()])-1):
          lowedge = requiredbinning[era][njet+"_"+ptll][var.lower()][i]
          highedge = requiredbinning[era][njet+"_"+ptll][var.lower()][i+1]
          assert f"{lowedge}_{highedge}" in finalweight[era][njet+"_"+ptll][var.lower()]
          varbin["value"]["content"].append(
            { "nodetype": "formula",
              "expression": finalweight[era][njet+"_"+ptll][var.lower()][f"{lowedge}_{highedge}"],
              "parser": "TFormula",
              "variables": ["val"]
            })
        ptllbin["content"].append(varbin)
      jetbin["content"].append(ptllbin)
    orderbin["value"]["content"].append(jetbin)
  recera["value"]["content"].append(orderbin)
  byeraREC["data"]["content"].append(orderbin)
 correctionlibs["Full"][actera]["corrections"].append(byeraREC)
 correctionlibs["REC"][actera]["corrections"].append(byeraREC)
 REC["data"]["content"].append(recera)

correctionlibs["Full"]["AllEras"]["corrections"].append(REC)
correctionlibs["REC"]["AllEras"]["corrections"].append(REC)

################################################## RECOIL CORR UNCERTAINTY

alleras = list(dict.fromkeys([f'{era.split("_")[0]}' for era in recoildata['ee']]))
allorders = sorted(list(dict.fromkeys([f'{era.split("_")[1]}' for era in recoildata['ee']])))
REC = {}
REC["name"] = "Recoil_correction_Uncertainty"
REC["description"] = "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_Uncertainty' gives uncertainties which are to be applies in addition to the nominal corrections. Note that a different input is required. The uncertainties were derived based on corrections applied with the QuantileMapHist method."
REC["version"] = 1
REC["inputs"] = []
REC["inputs"].append({"name": "era", "type": "string", "description": f"Era: {', '.join(alleras)}"})
REC["inputs"].append({"name": "order", "type": "string", "description": f"Order of samples: {', '.join(allorders)}"})
REC["inputs"].append({"name": "njet", "type": "real", "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"})
REC["inputs"].append({"name": "ptll", "type": "real", "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"})
REC["inputs"].append({"name": "var", "type": "string", "description": "The variable name you are giving the input for: 'Hpara', 'Hperp' (string). Note that this is a DIFFERENT vairable from Upara/Uperp! The output will be for the same kind of variable."})
REC["inputs"].append({"name": "val", "type": "real", "description": "Input value of either Hpara or Hperp"})
REC["inputs"].append({"name": "syst", "type": "string", "description": "Systematic variations on Response and Resolution: 'RespUp', 'RespDown', 'ResolUp', 'ResolDown'"})
REC["output"] = {"name": "outval", "type": "real", "description": "Output value"}
REC["data"] = {}
REC["data"]["nodetype"] = "category"
REC["data"]["input"] = "era"
REC["data"]["content"] = []
plain = copy.deepcopy(REC)

for actera in alleras:
 recera = {}
 recera["key"] = actera
 recera["value"] = {}
 recera["value"]["nodetype"] = "category"
 recera["value"]["input"] = "order"
 recera["value"]["content"] = []
 byeraREC = copy.deepcopy(plain)
 del byeraREC["inputs"][0]
 byeraREC["data"]["input"] = "order"
 for order in allorders:
  era = f'{actera}_{order}'
  if era not in recoildata["ee"]: continue
  orderbin = {}
  orderbin["key"] = order
  orderbin["value"] = {}
  orderbin["value"]["nodetype"] = "binning"
  orderbin["value"]["flow"] = "error"
  orderbin["value"]["input"] = "njet"
  orderbin["value"]["edges"] = [-0.5, 0.5, 1.5, float("inf")]
  orderbin["value"]["content"] = []
  for njet in ["0j", "1j", "2j"]:
    jetbin = {}
    jetbin["nodetype"] = "binning"
    jetbin["flow"] = "error"
    jetbin["input"] = "ptll"
    jetbin["edges"] = [0.0, 5.0, 15.0, 25.0, 40.0, 75.0, float("inf")]
    jetbin["content"] = []
    for ptll in ["0-5", "5-15", "15-25", "25-40", "40-75", "75-Inf"]:
      ptllbin = {}
      ptllbin["nodetype"] = "category"
      ptllbin["input"] = "var"
      ptllbin["content"] = []
      for var in ["Hpara", "Hperp"]:
        varbin = {}
        varbin["key"] = var
        varbin["value"] = {}
        varbin["value"]["nodetype"] = "category"
        varbin["value"]["input"] = "syst"
        varbin["value"]["content"] = []
        for sys in ["RespUp", "RespDown", "ResolUp", "ResolDown"]:
          sysbin = {}
          IsResp = ("Resp" in sys)
          IsUp = 1 if "Up" in sys else -1
          sysbin["key"] = sys
          sysbin["value"] = {}
          sysbin["value"]["nodetype"] = "formula"
          if IsResp and var=="Hpara":
            sysbin["value"]["expression"] = f'x + ({IsUp}) * {unc[era]["Uncertainty"][njet]["response"]} * {unc[era]["Response"][njet+"_"+ptll]} * y'
          elif IsResp and var=="Hperp":
            sysbin["value"]["expression"] = f'x'
          elif not IsResp and var=="Hpara":
            sysbin["value"]["expression"] = f'{unc[era]["Response"][njet+"_"+ptll]} * y + (1 + ({IsUp}) * {unc[era]["Uncertainty"][njet]["resolution"]}) * (x - {unc[era]["Response"][njet+"_"+ptll]} * y)'
          elif not IsResp and var=="Hperp":
            sysbin["value"]["expression"] = f'x * (1 + ({IsUp}) * {unc[era]["Uncertainty"][njet]["resolution"]})'
          sysbin["value"]["parser"] = "TFormula"
          sysbin["value"]["variables"] = ["val", "ptll"]
          varbin["value"]["content"].append(sysbin)
        ptllbin["content"].append(varbin)
      jetbin["content"].append(ptllbin)
    orderbin["value"]["content"].append(jetbin)
  recera["value"]["content"].append(orderbin)
  byeraREC["data"]["content"].append(orderbin)
 correctionlibs["Full"][actera]["corrections"].append(byeraREC)
 correctionlibs["REC"][actera]["corrections"].append(byeraREC)
 REC["data"]["content"].append(recera)

correctionlibs["Full"]["AllEras"]["corrections"].append(REC)
correctionlibs["REC"]["AllEras"]["corrections"].append(REC)

### Output

for era in correctionlibs["Full"]:
  erastring = "" if era=="AllEras" else "_"+era
  SaveCorrLib(correctionlibs["Full"][era], "DY_pTll_recoil_corrections"+erastring+".json")
for era in correctionlibs["DY"]:
  erastring = "" if era=="AllEras" else "_"+era
  SaveCorrLib(correctionlibs["DY"][era], "DY_pTll_weights"+erastring+".json")
for era in correctionlibs["REC"]:
  erastring = "" if era=="AllEras" else "_"+era
  SaveCorrLib(correctionlibs["REC"][era], "Recoil_corrections"+erastring+".json")

print("Done!")
exit()





