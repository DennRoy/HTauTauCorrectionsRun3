# HTauTauCorrectionsRun3

Scripts that I used to obtain corrections on DY MC over pTll, as well as recoil corrections. Another script combines these results and puts them in a file that can be read in correctionlib format.

# DYpTllCorr

Setup used to obtain pTll weights. The dictionary defined in the very beginning of the script reads histograms from rootfiles, which were produced in a previous step with a different framework. These histograms are distributions of events passing selection criteria (like di-muon final state) for different processes (Data, DY, other backgrounds). The script does all the work for each set of inputs separately.

# RecoilCorr

Setup used to obtain recoil corrections and uncertainties. Similar to the pTll script, "RecoilCorr.py" orks with loaded histograms from rootfiles produced beforehand. The inputs defined with the appendix "_appl" are those with have previously obtained recoil corrections already applied, and those are needed to obtain uncertainties.

# TauCorrLib

Contains the script which takes all outputs produced by `DYpTllCorr` and `RecoilCorr` as input, to produce a correctionlib file. Includes some complex code on inverting double gaussian functions that have been split into many linear functions.
