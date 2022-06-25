#!/usr/bin/env python
import ROOT
import imp
import json
from array import array
import math
import re
wsptools = imp.load_source('wsptools', 'workspaceTools.py')


def GetFromTFile(str):
    f = ROOT.TFile(str.split(':')[0])
    obj = f.Get(str.split(':')[1]).Clone()
    f.Close()
    return obj


# Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.TH1.AddDirectory(0)
ROOT.gROOT.LoadMacro("CrystalBallEfficiency.cxx+")

w = ROOT.RooWorkspace('w')

########################################
### Muon trk efficiency from MuonPOG ###
########################################
loc = 'inputs/2016UL_postVFP/MuonPOG'

histsToWrap = [
    (loc+'/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root:NUM_TrackerMuons_DEN_genTracks',
      'm_trk_eff_lowpt'),
    (loc+'/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt.root:NUM_TrackerMuons_DEN_genTracks',
     'm_trk_eff_medpt')
]
for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['m_eta', 'm_pt'],
                          GetFromTFile(task[0]), name=task[1])

wsptools.MakeBinnedCategoryFuncMap(w, 'm_pt', [3., 15., 120.],
                                   'm_trk_ratio', ['m_trk_eff_lowpt', 'm_trk_eff_medpt'])
                                   
###############################################
### Electron reco efficiency from EGammaPOG ###
###############################################
loc = 'inputs/2016UL_postVFP/EGammaPOG'

histsToWrap = [
    (loc+'/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D',
     'e_trk_ST20_ratio'),
    (loc+'/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D',
     'e_trk_GT20_ratio')
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['e_eta', 'e_pt'],
                          GetFromTFile(task[0]), name=task[1])

wsptools.MakeBinnedCategoryFuncMap(w, 'e_pt', [10., 20., 500.],
                                   'e_trk_ratio', ['e_trk_ST20_ratio', 'e_trk_GT20_ratio'])

#######################################
### deepTau Trigger SFs from TauPOG ###
#######################################

loc = 'inputs/2016UL_postVFP/TauPOGTrigger/'
tau_trg_file = ROOT.TFile(loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root')
#tau_id_wps=['VVVLoose','VVLoose','VLoose','Loose','Medium','Tight']
tau_id_wps=['Medium']#,'Tight']

for wp in tau_id_wps:
  for dm in ['0','1','10','11']:
    histsToWrap = [
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:data_ditau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_ditau_dm%s_data' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:mc_ditau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_ditau_dm%s_mc' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:sf_ditau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_ditau_dm%s_ratio' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:data_mutau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_mutau_dm%s_data' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:mc_mutau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_mutau_dm%s_mc' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:sf_mutau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_mutau_dm%s_ratio' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:data_etau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_etau_dm%s_data' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:mc_etau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_etau_dm%s_mc' % (wp.lower(),dm)),
      (loc+'2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root:sf_etau_%s_dm%s_fitted' % (wp,dm),  't_trg_pog_deeptau_%s_etau_dm%s_ratio' % (wp.lower(),dm)),
    ]

    for task in histsToWrap:
        wsptools.SafeWrapHist(w, ['t_pt'],
                              GetFromTFile(task[0]), name=task[1])

        hist = GetFromTFile(task[0])
        uncert_hists = wsptools.UncertsFromHist(hist)
        wsptools.SafeWrapHist(w, ['t_pt'], uncert_hists[0], name=task[1]+'_up')
        wsptools.SafeWrapHist(w, ['t_pt'], uncert_hists[1], name=task[1]+'_down')

        if 'ditau' in task[1]:
          wsptools.SafeWrapHist(w, ['t_pt_2'],
                                GetFromTFile(task[0]), name=task[1]+'_2')
  
          hist = GetFromTFile(task[0])
          uncert_hists = wsptools.UncertsFromHist(hist)
          wsptools.SafeWrapHist(w, ['t_pt_2'], uncert_hists[0], name=task[1]+'_up_2')
          wsptools.SafeWrapHist(w, ['t_pt_2'], uncert_hists[1], name=task[1]+'_down_2') 

  wp_lower = wp.lower()
  for i in ['data','mc','ratio']:
    for j in ['ditau','mutau', 'etau']:
      taus=['']
      if j == 'ditau': taus = ['', '_2']
      for t in taus:
        w.factory('expr::t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s%(t)s("(@0==0)*@1 + (@0==1||@0==2)*@2 + (@0==10||@0==5)*@3 + (@0==11||@0==6)*@4", t_dm%(t)s[0], t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm0_%(i)s%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm1_%(i)s%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm10_%(i)s%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm11_%(i)s%(t)s)' % vars())

        w.factory('expr::t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_up%(t)s("@5 + ((@0==0)*@1 + (@0==1||@0==2)*@2 + (@0==10||@0==5)*@3 + (@0==11||@0==6)*@4)", t_dm%(t)s[0], t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm0_%(i)s_up%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm1_%(i)s_up%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm10_%(i)s_up%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm11_%(i)s_up%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s%(t)s)' % vars())

        w.factory('expr::t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_down%(t)s("@5 - ((@0==0)*@1 + (@0==1||@0==2)*@2 + (@0==10||@0==5)*@3 + (@0==11||@0==6)*@4)", t_dm%(t)s[0], t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm0_%(i)s_down%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm1_%(i)s_down%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm10_%(i)s_down%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_dm11_%(i)s_down%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s%(t)s)' % vars())

        for dm in ['0','1','10','11']:
          w.factory('expr::t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_dm%(dm)s_down%(t)s("(@0==%(dm)s)*@1 + (@0!=%(dm)s)*@2",t_dm%(t)s[0], t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_down%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s%(t)s)' % vars())
          w.factory('expr::t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_dm%(dm)s_up%(t)s("(@0==%(dm)s)*@1 + (@0!=%(dm)s)*@2",t_dm%(t)s[0], t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s_up%(t)s, t_trg_pog_deeptau_%(wp_lower)s_%(j)s_%(i)s%(t)s)' % vars())

##################################
### deepTau ID SFs from TauPOG ###
##################################
w.factory('expr::t_dm_bounded("(@0<2)*@0 +(@0==2)*1 + (@0>2&&@0<11)*10 + (@0>10)*11" ,t_dm[0])')

# dm binned SFs
loc = 'inputs/2016UL_postVFP/TauPOGID'

histsToWrap = [
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:VVVLoose', 't_deeptauid_dm_vvvloose'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:VVLoose', 't_deeptauid_dm_vvloose'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:VLoose', 't_deeptauid_dm_vloose'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:Loose', 't_deeptauid_dm_loose'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:Medium', 't_deeptauid_dm_medium'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:Tight', 't_deeptauid_dm_tight'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:VTight', 't_deeptauid_dm_vtight'),
  (loc+'/TauID_SF_dm_DeepTau2017v2p1VSjet_UL2016_postVFP.root:VVTight', 't_deeptauid_dm_vvtight'),
]

for task in histsToWrap:
  wsptools.SafeWrapHist(w, ['t_dm_bounded'], GetFromTFile(task[0]), name=task[1])
  uncert_hists = wsptools.UncertsFromHist(GetFromTFile(task[0]))
  wsptools.SafeWrapHist(w, ['t_dm_bounded'], uncert_hists[0], name=task[1]+'_abs_up')
  wsptools.SafeWrapHist(w, ['t_dm_bounded'], uncert_hists[1], name=task[1]+'_abs_down')
  w.factory('expr::%s_up("@1+@0",%s_abs_up,%s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_down("@1-@0",%s_abs_down,%s)' % (task[1],task[1],task[1]))

  w.factory('expr::%s_dm0_up("(@0==0)*@1 + (@0!=0)*@2 ", t_dm[0], %s_up, %s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_dm0_down("(@0==0)*@1 + (@0!=0)*@2 ", t_dm[0], %s_down, %s)' % (task[1],task[1],task[1]))

  w.factory('expr::%s_dm1_up("(@0==1)*@1 + (@0!=1)*@2 ", t_dm[0], %s_up, %s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_dm1_down("(@0==1)*@1 + (@0!=1)*@2 ", t_dm[0], %s_down, %s)' % (task[1],task[1],task[1]))

  w.factory('expr::%s_dm10_up("(@0==10)*@1 + (@0!=10)*@2 ", t_dm[0], %s_up, %s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_dm10_down("(@0==10)*@1 + (@0!=10)*@2 ", t_dm[0], %s_down, %s)' % (task[1],task[1],task[1]))

  w.factory('expr::%s_dm11_up("(@0==11)*@1 + (@0!=11)*@2 ", t_dm[0], %s_up, %s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_dm11_down("(@0==11)*@1 + (@0!=11)*@2 ", t_dm[0], %s_down, %s)' % (task[1],task[1],task[1]))

# pT dependent SFs (TauID)
sf_funcs = {}
tauid_pt_file = ROOT.TFile(loc+'/TauID_SF_pt_DeepTau2017v2p1VSjet_UL2016_postVFP.root')
for i in ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
  for j in ['cent', 'up', 'down']:
    fname = '%s_%s' % (i,j)
    fit = tauid_pt_file.Get(fname)
    outname = i.lower()
    if j != 'cent': outname+='_%s' % j
    sf_funcs[outname] = fit.GetTitle()

for x in sf_funcs:
  func = re.sub('x','@0',sf_funcs[x])
  w.factory('expr::t_deeptauid_pt_%s("%s",t_pt[0])' % (x, func))
# pt_bins [20,25,30,35,40,infty]
for i in ['vvvloose', 'vvloose', 'vloose', 'loose', 'medium', 'tight', 'vtight', 'vvtight']:
  w.factory('expr::t_deeptauid_pt_%(i)s_bin1_up("(@0>20&&@0<=25)*@1 + ((@0>20&&@0<=25)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin1_down("(@0>20&&@0<=25)*@1 + ((@0>20&&@0<=25)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_pt_%(i)s_bin2_up("(@0>25&&@0<=30)*@1 + ((@0>25&&@0<=30)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin2_down("(@0>25&&@0<=30)*@1 + ((@0>25&&@0<=30)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_pt_%(i)s_bin3_up("(@0>30&&@0<=35)*@1 + ((@0>30&&@0<=35)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin3_down("(@0>30&&@0<=35)*@1 + ((@0>30&&@0<=35)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_pt_%(i)s_bin4_up("(@0>35&&@0<=40)*@1 + ((@0>35&&@0<=40)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin4_down("(@0>35&&@0<=40)*@1 + ((@0>35&&@0<=40)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_pt_%(i)s_bin5_up("(@0>40&&@0<=500)*@1 + ((@0>40&&@0<=500)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin5_down("(@0>40&&@0<=500)*@1 + ((@0>40&&@0<=500)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  #w.factory('expr::t_deeptauid_pt_%(i)s_bin5_up("(@0>40)*@1 + ((@0>40)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  #w.factory('expr::t_deeptauid_pt_%(i)s_bin5_down("(@0>40)*@1 + ((@0>40)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_pt_%(i)s_bin6_up("(@0>500)*@1 + ((@0>500)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_up, t_deeptauid_pt_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_pt_%(i)s_bin6_down("(@0>500)*@1 + ((@0>500)==0)*@2",t_pt[0], t_deeptauid_pt_%(i)s_down, t_deeptauid_pt_%(i)s)' % vars())

  w.factory('expr::t_deeptauid_dm_%(i)s_bin6_up("((@0>500)*@1/@2 + ((@0>500)==0))*@3",t_pt[0], t_deeptauid_pt_%(i)s_bin6_up, t_deeptauid_pt_%(i)s, t_deeptauid_dm_%(i)s)' % vars())
  w.factory('expr::t_deeptauid_dm_%(i)s_bin6_down("((@0>500)*@1/@2 + ((@0>500)==0))*@3",t_pt[0], t_deeptauid_pt_%(i)s_bin6_down, t_deeptauid_pt_%(i)s, t_deeptauid_dm_%(i)s)' % vars())

# l->tau fake scale factors - Updated for Ultra Legacy Samples
loc='inputs/2016UL_postVFP/TauPOGID/'

histsToWrap = [
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:VVLoose', 't_id_vs_e_eta_vvloose'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:VLoose', 't_id_vs_e_eta_vloose'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:Loose',  't_id_vs_e_eta_loose'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:Medium', 't_id_vs_e_eta_medium'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:Tight',  't_id_vs_e_eta_tight'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:VTight', 't_id_vs_e_eta_vtight'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSe_UL2016_postVFP.root:VVTight', 't_id_vs_e_eta_vvtight'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root:VLoose', 't_id_vs_mu_eta_vloose'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root:Loose',  't_id_vs_mu_eta_loose'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root:Medium', 't_id_vs_mu_eta_medium'),
  (loc+'/TauID_SF_eta_DeepTau2017v2p1VSmu_UL2016_postVFP.root:Tight',  't_id_vs_mu_eta_tight'),
]

w.factory('expr::t_eta_bounded("min(2.3,TMath::Abs(@0))" ,t_eta[0])')

for task in histsToWrap:
  wsptools.SafeWrapHist(w, ['t_eta_bounded'], GetFromTFile(task[0]), name=task[1])
  uncert_hists = wsptools.UncertsFromHist(GetFromTFile(task[0]))
  wsptools.SafeWrapHist(w, ['t_eta_bounded'], uncert_hists[0], name=task[1]+'_abs_up')
  wsptools.SafeWrapHist(w, ['t_eta_bounded'], uncert_hists[1], name=task[1]+'_abs_down')
  w.factory('expr::%s_up("@1+@0",%s_abs_up,%s)' % (task[1],task[1],task[1]))
  w.factory('expr::%s_down("@1-@0",%s_abs_down,%s)' % (task[1],task[1],task[1]))


###############################################
### IC electron id, iso, trigger SFs for MC ###
###############################################

loc = 'inputs/2016UL_postVFP/ICSF/tpzee/'

histsToWrap = [

    (loc+'singleElec/electron_SFs.root:data_trg_eff', 'e_trg_ic_data'),
    (loc+'singleElec/electron_SFs.root:ZLL_trg_eff', 'e_trg_ic_mc'),

    (loc+'singleElec/electron_SFs.root:data_iso_eff', 'e_iso_ic_data'),
    (loc+'singleElec/electron_SFs.root:ZLL_iso_eff', 'e_iso_ic_mc'),

    (loc+'singleElec/electron_SFs.root:data_id_eff', 'e_id_ic_data'),
    (loc+'singleElec/electron_SFs.root:ZLL_id_eff', 'e_id_ic_mc'),
    
    (loc+'emLow/electron_SFs.root:data_trg_eff', 'e_trg_12_ic_data'),
    (loc+'emLow/electron_SFs.root:ZLL_trg_eff', 'e_trg_12_ic_mc'),
    
    (loc+'emHigh/electron_SFs.root:data_trg_eff', 'e_trg_23_ic_data'),
    (loc+'emHigh/electron_SFs.root:ZLL_trg_eff', 'e_trg_23_ic_mc'),
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])'],
                          GetFromTFile(task[0]), name=task[1])

# temporarily take isolated SFs for all until anti iso ones are included
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_trg_binned_ic_data', ['e_trg_ic_data', 'e_trg_ic_data', 'e_trg_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_trg_binned_ic_mc', ['e_trg_ic_mc', 'e_trg_ic_mc', 'e_trg_ic_mc'])
                                   
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_iso_binned_ic_data', ['e_iso_ic_data', 'e_iso_ic_data', 'e_iso_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_iso_binned_ic_mc', ['e_iso_ic_mc', 'e_iso_ic_mc', 'e_iso_ic_mc'])
                                   
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.50],
                                   'e_trg_12_binned_ic_data', ['e_trg_12_ic_data', 'e_trg_12_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.50],
                                   'e_trg_12_binned_ic_mc', ['e_trg_12_ic_mc', 'e_trg_12_ic_mc'])

wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.50],
                                   'e_trg_23_binned_ic_data', ['e_trg_23_ic_data', 'e_trg_23_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.50],
                                   'e_trg_23_binned_ic_mc', ['e_trg_23_ic_mc', 'e_trg_23_ic_mc'])
                                   
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_trg_24_binned_ic_data', ['e_trg_24_ic_data', 'e_trg_24_ic_data', 'e_trg_24_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'e_iso', [0., 0.15, 0.25, 0.50],
                                   'e_trg_24_binned_ic_mc', ['e_trg_24_ic_mc', 'e_trg_24_ic_mc', 'e_trg_24_ic_mc'])                                   
                                   
                                                              
w.factory('expr::e_idiso_ic_data("@0*@1", e_iso_ic_data, e_id_ic_data)' % vars())
w.factory('expr::e_idiso_ic_mc("@0*@1", e_iso_ic_mc, e_id_ic_mc)' % vars())

w.factory('expr::e_idiso_binned_ic_data("@0*@1", e_iso_binned_ic_data, e_id_ic_data)' % vars())
w.factory('expr::e_idiso_binned_ic_mc("@0*@1", e_iso_binned_ic_mc, e_id_ic_mc)' % vars())

for i in ['trg', 'id', 'iso', 'idiso','trg_12','trg_23']:
  w.factory('expr::e_%(i)s_ic_ratio("@0/@1", e_%(i)s_ic_data, e_%(i)s_ic_mc)' % vars())
  w.factory('expr::e_%(i)s_binned_ic_ratio("@0/@1", e_%(i)s_binned_ic_data, e_%(i)s_binned_ic_mc)' % vars())
# ask about id

###########################################
### IC muon id, iso, trigger SFs for MC ###
###########################################

loc = 'inputs/2016UL_postVFP/ICSF/tpzmm/'

histsToWrap = [
    (loc+'singleMu/muon_SFs.root:data_trg_eff', 'm_trg_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_trg_eff', 'm_trg_ic_mc'),

    (loc+'singleMu/muon_SFs.root:data_iso_eff', 'm_iso_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_iso_eff', 'm_iso_ic_mc'),

    (loc+'singleMu/muon_SFs.root:data_id_eff', 'm_id_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_id_eff', 'm_id_ic_mc'),
        
    #mt
    (loc+'mt/muon_SFs.root:data_trg_eff', 'm_trg_19_ic_data'),
    (loc+'mt/muon_SFs.root:ZLL_trg_eff', 'm_trg_19_ic_mc'),
        
    #em high and low legs
    (loc+'emLow/muon_SFs.root:data_trg_eff', 'm_trg_8_ic_data'),
    (loc+'emLow/muon_SFs.root:ZLL_trg_eff', 'm_trg_8_ic_mc'),
    
    (loc+'emLow/muon_SFs.root:data_iso_eff', 'm_looseiso_ic_data'),
    (loc+'emLow/muon_SFs.root:ZLL_iso_eff', 'm_looseiso_ic_mc'),
    
    (loc+'emHigh/muon_SFs.root:data_trg_eff', 'm_trg_23_ic_data'),
    (loc+'emHigh/muon_SFs.root:ZLL_trg_eff', 'm_trg_23_ic_mc'),

    (loc+'emLow/muon_SFs.root:data_id_eff', 'm_id_ic_data'),
    (loc+'emLow/muon_SFs.root:ZLL_id_eff', 'm_id_ic_mc'),
        
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)",m_eta[0])'],
                          GetFromTFile(task[0]), name=task[1])


wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_binned_ic_data', ['m_trg_ic_data', 'm_trg_ic_data', 'm_trg_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_binned_ic_mc', ['m_trg_ic_mc', 'm_trg_ic_mc', 'm_trg_ic_mc'])

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_iso_binned_ic_data', ['m_iso_ic_data', 'm_iso_ic_data', 'm_iso_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_iso_binned_ic_mc', ['m_iso_ic_mc', 'm_iso_ic_mc', 'm_iso_ic_mc'])                                   

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_trg_8_binned_ic_data', ['m_trg_8_ic_data', 'm_trg_8_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_trg_8_binned_ic_mc', ['m_trg_8_ic_mc', 'm_trg_8_ic_mc'])

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_looseiso_binned_ic_data', ['m_looseiso_ic_data', 'm_looseiso_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_looseiso_binned_ic_mc', ['m_looseiso_ic_mc', 'm_looseiso_ic_mc'])

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_trg_23_binned_ic_data', ['m_trg_23_ic_data', 'm_trg_23_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.2, 0.50],
                                   'm_trg_23_binned_ic_mc', ['m_trg_23_ic_mc', 'm_trg_23_ic_mc'])
                                
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_19_binned_ic_data', ['m_trg_19_ic_data', 'm_trg_19_ic_data', 'm_trg_19_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_19_binned_ic_mc', ['m_trg_19_ic_mc', 'm_trg_19_ic_mc', 'm_trg_19_ic_mc'])



w.factory('expr::m_idiso_ic_data("@0*@1", m_iso_ic_data, m_id_ic_data)' % vars())
w.factory('expr::m_idiso_ic_mc("@0*@1", m_iso_ic_mc, m_id_ic_mc)' % vars())

w.factory('expr::m_idiso_binned_ic_data("@0*@1", m_iso_binned_ic_data, m_id_ic_data)' % vars())
w.factory('expr::m_idiso_binned_ic_mc("@0*@1", m_iso_binned_ic_mc, m_id_ic_mc)' % vars())

w.factory('expr::m_idlooseiso_ic_data("@0*@1", m_looseiso_ic_data, m_id_ic_data)' % vars())
w.factory('expr::m_idlooseiso_ic_mc("@0*@1", m_looseiso_ic_mc, m_id_ic_mc)' % vars())

w.factory('expr::m_idlooseiso_binned_ic_data("@0*@1", m_looseiso_binned_ic_data, m_id_ic_data)' % vars())
w.factory('expr::m_idlooseiso_binned_ic_mc("@0*@1", m_looseiso_binned_ic_mc, m_id_ic_mc)' % vars())

for i in ['trg', 'trg_19', 'trg_8', 'trg_23', 'id', 'iso', 'looseiso', 'idiso', 'idlooseiso']:
  w.factory('expr::m_%(i)s_ic_ratio("@0/@1", m_%(i)s_ic_data, m_%(i)s_ic_mc)' % vars())
  w.factory('expr::m_%(i)s_ic_embed_ratio("@0/@1", m_%(i)s_ic_data, m_%(i)s_ic_embed)' % vars())
  w.factory('expr::m_%(i)s_binned_ic_ratio("@0/@1", m_%(i)s_binned_ic_data, m_%(i)s_binned_ic_mc)' % vars())
  w.factory('expr::m_%(i)s_binned_ic_embed_ratio("@0/@1", m_%(i)s_binned_ic_data, m_%(i)s_binned_ic_embed)' % vars())


histsToWrap = [
    (loc+'dimuLow/muon_SFs.root:data_id_eff', 'm_sel_id_ic_1_data'),
    (loc+'dimuLow/muon_SFs.root:data_trg_eff', 'm_sel_trg_8_ic_1_data'),
    (loc+'dimuHigh/muon_SFs.root:data_trg_eff', 'm_sel_trg_17_ic_1_data'),
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['gt1_pt', 'expr::gt1_abs_eta("TMath::Abs(@0)",gt1_eta[0])'],
                          GetFromTFile(task[0]), name=task[1])

histsToWrap = [
    (loc+'dimuLow/muon_SFs.root:data_id_eff', 'm_sel_id_ic_2_data'),
    (loc+'dimuLow/muon_SFs.root:data_trg_eff', 'm_sel_trg_8_ic_2_data'),
    (loc+'dimuHigh/muon_SFs.root:data_trg_eff', 'm_sel_trg_17_ic_2_data'),
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['gt2_pt', 'expr::gt2_abs_eta("TMath::Abs(@0)",gt2_eta[0])'],
                          GetFromTFile(task[0]), name=task[1])

w.factory('expr::m_sel_trg_ic_data("0.9454*(@0*@3+@1*@2-@1*@3)", m_sel_trg_8_ic_1_data, m_sel_trg_17_ic_1_data, m_sel_trg_8_ic_2_data, m_sel_trg_17_ic_2_data)')
w.factory('expr::m_sel_trg_ic_ratio("min(1./@0,20)", m_sel_trg_ic_data)')

wsptools.SafeWrapHist(w, ['gt_pt', 'expr::gt_abs_eta("TMath::Abs(@0)",gt_eta[0])'],
                          GetFromTFile(loc+'dimuLow/muon_SFs.root:data_id_eff'), 'm_sel_id_ic_data')

w.factory('expr::m_sel_id_ic_ratio("min(1./@0,20)", m_sel_id_ic_data)')


w.importClassCode('CrystalBallEfficiency')
w.Print()
w.writeToFile('output/htt_scalefactors_UL_2016postVFP.root')
w.Delete()
