#!/usr/bin/env python
import ROOT
import imp
import json
from array import array
import numpy as np
import math
import re
import os

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

###########################################
### IC muon id, iso, trigger SFs for MC ###
###########################################
loc = 'inputs/2022_postEE/ICSF/tpzmm/'

histsToWrap = [

    (loc+'singleMu/muon_SFs.root:data_id_eff', 'm_id_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_id_eff', 'm_id_ic_mc'),
    #(loc+'singleMu/muon_SFs.root:embed_id_eff', 'm_id_ic_embed'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_EMB_id', 'm_id_ic_embed_ratio'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_id', 'm_id_ic_ratio'),

    (loc+'singleMu/muon_SFs.root:data_iso_eff', 'm_iso_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_iso_eff', 'm_iso_ic_mc'),
    #(loc+'singleMu/muon_SFs.root:embed_iso_eff', 'm_iso_ic_embed'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_EMB_iso', 'm_iso_ic_embed_ratio'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_iso', 'm_iso_ic_ratio'),

    (loc+'singleMu/muon_SFs.root:data_trg_eff', 'm_trg_ic_data'),
    (loc+'singleMu/muon_SFs.root:ZLL_trg_eff', 'm_trg_ic_mc'),
    #(loc+'singleMu/muon_SFs.root:embed_trg_eff', 'm_trg_ic_embed'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_EMB_trg', 'm_trg_ic_embed_ratio'),
    #(loc+'singleMu/muon_SFs.root:ScaleFactor_trg', 'm_trg_ic_ratio'),
    
#    (loc+'emLow/muon_SFs.root:data_trg_eff', 'm_trg_8_ic_data'),
#    (loc+'emLow/muon_SFs.root:ZLL_trg_eff', 'm_trg_8_ic_mc'),
    
#    (loc+'emLow/muon_SFs.root:data_iso_eff', 'm_looseiso_ic_data'),
#    (loc+'emLow/muon_SFs.root:ZLL_iso_eff', 'm_looseiso_ic_mc'),
    
#    (loc+'emHigh/muon_SFs.root:data_trg_eff', 'm_trg_23_ic_data'),
#    (loc+'emHigh/muon_SFs.root:ZLL_trg_eff', 'm_trg_23_ic_mc'),
]
                          
for task in histsToWrap:
   wsptools.SafeWrapHist(w, ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)",m_eta[0])'], GetFromTFile(task[0]), name=task[1])
   uncert_hists = wsptools.UncertsFromHist2D(GetFromTFile(task[0]))
   wsptools.SafeWrapHist(w, ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)",m_eta[0])'], uncert_hists[0], name=task[1]+'_abs_up')
   wsptools.SafeWrapHist(w, ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)",m_eta[0])'], uncert_hists[1], name=task[1]+'_abs_down')
   w.factory('expr::%s_up("@1+@0",%s_abs_up,%s)' % (task[1],task[1],task[1]))
   w.factory('expr::%s_down("@1-@0",%s_abs_down,%s)' % (task[1],task[1],task[1]))

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_binned_ic_data', ['m_trg_ic_data', 'm_trg_ic_data', 'm_trg_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_trg_binned_ic_mc', ['m_trg_ic_mc', 'm_trg_ic_mc', 'm_trg_ic_mc'])
#wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
#                                   'm_trg_binned_ic_embed', ['m_trg_ic_embed', 'm_trg_ic_embed', 'm_trg_ic_embed'])

#wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
#                                   'm_trg_20_binned_ic_data', ['m_trg_20_ic_data', 'm_trg_20_ic_data', 'm_trg_20_ic_data'])
#wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
#                                   'm_trg_20_binned_ic_mc', ['m_trg_20_ic_mc', 'm_trg_20_ic_mc', 'm_trg_20_ic_mc'])
#wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
#                                   'm_trg_20_binned_ic_embed', ['m_trg_20_ic_embed', 'm_trg_20_ic_embed', 'm_trg_20_ic_embed'])

wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_iso_binned_ic_data', ['m_iso_ic_data', 'm_iso_ic_data', 'm_iso_ic_data'])
wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
                                   'm_iso_binned_ic_mc', ['m_iso_ic_mc', 'm_iso_ic_mc', 'm_iso_ic_mc'])
#wsptools.MakeBinnedCategoryFuncMap(w, 'm_iso', [0., 0.15, 0.25, 0.50],
#                                   'm_iso_binned_ic_embed', ['m_iso_ic_embed', 'm_iso_ic_embed', 'm_iso_ic_embed'])

w.factory('expr::m_idiso_ic_data("@0*@1", m_id_ic_data, m_iso_ic_data)')
w.factory('expr::m_idiso_ic_mc("@0*@1", m_id_ic_mc, m_iso_ic_mc)')
#w.factory('expr::m_idiso_ic_embed("@0*@1", m_id_ic_embed, m_iso_ic_embed)')
#w.factory('expr::m_idiso_ic_ratio("@0*@1", m_id_ic_ratio, m_iso_ic_ratio)')
#w.factory('expr::m_idiso_ic_embed_ratio("@0*@1", m_id_ic_embed_ratio, m_iso_ic_embed_ratio)')
# w.factory('expr::m_idlooseiso_ic_data("@0*@1", m_looseiso_ic_data, m_id_ic_data)' % vars())
# w.factory('expr::m_idlooseiso_ic_mc("@0*@1", m_looseiso_ic_mc, m_id_ic_mc)' % vars())

w.factory('expr::m_idiso_binned_ic_data("@0*@1", m_iso_binned_ic_data, m_id_ic_data)' % vars())
w.factory('expr::m_idiso_binned_ic_mc("@0*@1", m_iso_binned_ic_mc, m_id_ic_mc)' % vars())
#w.factory('expr::m_idiso_binned_ic_embed("@0*@1", m_iso_binned_ic_embed, m_id_ic_embed)' % vars())
# w.factory('expr::m_idlooseiso_binned_ic_data("@0*@1", m_looseiso_binned_ic_data, m_id_ic_data)' % vars())
# w.factory('expr::m_idlooseiso_binned_ic_mc("@0*@1", m_looseiso_binned_ic_mc, m_id_ic_mc)' % vars())
#w.factory('expr::m_idlooseiso_binned_ic_embed("@0*@1", m_looseiso_binned_ic_embed, m_id_ic_embed)' % vars())

for i in ['trg', 'iso', 'idiso']:#['trg', 'trg_20', 'trg_8', 'trg_23', 'id', 'iso', 'looseiso', 'idiso', 'idlooseiso']:
  w.factory('expr::m_%(i)s_ic_ratio("@0/@1", m_%(i)s_ic_data, m_%(i)s_ic_mc)' % vars())
#  w.factory('expr::m_%(i)s_ic_embed_ratio("@0/@1", m_%(i)s_ic_data, m_%(i)s_ic_embed)' % vars())
  w.factory('expr::m_%(i)s_binned_ic_ratio("@0/@1", m_%(i)s_binned_ic_data, m_%(i)s_binned_ic_mc)' % vars())
#  w.factory('expr::m_%(i)s_binned_ic_embed_ratio("@0/@1", m_%(i)s_binned_ic_data, m_%(i)s_binned_ic_embed)' % vars())
w.factory('expr::m_id_ic_ratio("@0/@1", m_id_ic_data, m_id_ic_mc)' % vars())

# for variation in ["up","down"]:
#       if variation == "up": var = "_up"
#       elif variation == "down": var = "_down"
#       w.factory('expr::m_idiso_ic_data{0}("TMath::Power(TMath::Power((@0/@1),2) + TMath::Power((@2/@3),2),0.5)*(@4),m_id_ic_data{0},m_id_ic_data,m_iso_ic_data{0},m_iso_ic_data,m_idiso_ic_data)'.format(var))
#       w.factory('expr::m_idiso_ic_mc{0}("TMath::Power(TMath::Power((@0/@1),2) + TMath::Power((@2/@3),2),0.5)*(@4),m_id_ic_mc{0},m_id_ic_mc,m_iso_ic_mc{0},m_iso_ic_mc,m_idiso_ic_mc)'.format(var))
#       #w.factory('expr::m_idiso_ic_embed{0}("TMath::Power(TMath::Power((@0/@1),2) + TMath::Power((@2/@3),2),0.5)*(@4),m_id_ic_embed{0},m_id_ic_embed,m_iso_ic_embed{0},m_iso_ic_embed,m_idiso_ic_embed)'.format(var))
#       w.factory('expr::m_idiso_ic_ratio{0}("TMath::Power(TMath::Power((@0/@1),2) + TMath::Power((@2/@3),2),0.5)*(@4),m_id_ic_ratio{0},m_id_ic_ratio,m_iso_ic_ratio{0},m_iso_ic_ratio,m_idiso_ic_ratio)'.format(var))
#       #w.factory('expr::m_idiso_ic_embed_ratio{0}("TMath::Power(TMath::Power((@0/@1),2) + TMath::Power((@2/@3),2),0.5)*(@4),m_id_ic_embed_ratio{0},m_id_ic_embed_ratio,m_iso_ic_embed_ratio{0},m_iso_ic_embed_ratio,m_idiso_ic_embed_ratio)'.format(var))
                                   

#zpt_reweighting LO & NLO
histsToWrap = [
    ('inputs/zpt/zpt_reweighting_LO_2022postEE.root:zptmass_histo', 'zptmass_weight_nom'),
    #('inputs/zpt/zpt_reweighting_NLO.root:zptmass_histo', 'zptmass_weight_nom_NLO')
]

for task in histsToWrap:
    wsptools.SafeWrapHist(w, ['z_gen_mass', 'z_gen_pt'],
                          GetFromTFile(task[0]), name=task[1])

############
### END  ###
############

w.importClassCode('CrystalBallEfficiency')

w.Print()

w.writeToFile('output/htt_scalefactors_2022postEE.root')

w.Delete()
