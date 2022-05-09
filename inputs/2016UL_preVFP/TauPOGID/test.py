import ROOT
import re 

sf_funcs = {}
tauid_pt_file = ROOT.TFile('TauID_SF_pt_DeepTau2017v2p1VSjet_UL2016_preVFP.root')
for i in ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
  for j in ['cent', 'up', 'down']:
    fname = '%s_%s' % (i,j)
    fit = tauid_pt_file.Get(fname)
    outname = i.lower()
    if j != 'cent': outname+='_%s' % j
    sf_funcs[outname] = fit.GetTitle()

for x in sf_funcs:
  func = re.sub('x','@0',sf_funcs[x])
  print()
  print(func)
