
import json
import re
import numpy as np
import ROOT
import sys

def read_json_sf(file):

  f = open(file)
  data = json.load(f)

  return data

def json_get_sf(identifier,sub_id, data):
  array = []
  eta_bins = []
  pt_bins = []
  full_arr = []
  for key in   data[identifier][sub_id].keys():
    if key != "binning":
      split_key =  re.findall(r'\[(.*?)\]', key)
      output = [x.split(',') for x in split_key]
      eta_bin_lower = float(output[0][0])
      eta_bin_upper = float(output[0][1])
      if eta_bin_lower not in eta_bins:
        eta_bins.append(eta_bin_lower)
      if eta_bin_upper not in eta_bins:
        eta_bins.append(eta_bin_upper)
      for sub_key in data[identifier][sub_id][key].keys():
        split_sub_key = re.findall(r'\[(.*?)\]', sub_key)
        sub_output = [x.split(',') for x in split_sub_key]
        pt_bin_lower = float(sub_output[0][0])
        pt_bin_upper = float(sub_output[0][1])
        temp_id = 0
        if pt_bin_lower not in pt_bins:
	  pt_bins.append(pt_bin_lower)
        if pt_bin_upper not in pt_bins:
	  pt_bins.append(pt_bin_upper)
        y =  data[identifier][sub_id][key][sub_key]["value"]
	if 'error' in data[identifier][sub_id][key][sub_key].keys():
          Ey =  data[identifier][sub_id][key][sub_key]["error"] 
          array.append([temp_id,eta_bin_lower,eta_bin_upper,temp_id,pt_bin_lower,pt_bin_upper,y,Ey])
        else: array.append([temp_id,eta_bin_lower,eta_bin_upper,temp_id,pt_bin_lower,pt_bin_upper,y])

  array = np.array(array)
  eta_bins = sorted(set(eta_bins))
  pt_bins = sorted(set(pt_bins))
  for i in range(len(array)):
    array[i][0] = eta_bins.index(array[i][1])
    array[i][3] = pt_bins.index(array[i][4])
  eta_bins = np.array(eta_bins)
  pt_bins = np.array(pt_bins)  
  return array,eta_bins,pt_bins

data = read_json_sf("%s.json" %sys.argv[1])
df,eta_bins,pt_bins = json_get_sf("%s" %sys.argv[2],"%s" %sys.argv[3],data)
# we have two custom binnings (eta_bins) and (pt_bins)
c = ROOT.TCanvas("c1","c1",800,1000)
h = ROOT.TH2F("h","%s"%sys.argv[2],len(eta_bins)-1,eta_bins,len(pt_bins)-1,pt_bins)

for i in df:
   
  h.SetBinContent(int(i[0]+1),int(i[3]+1),i[6])
  if sys.argv[4] == "do_errors":
     h.SetBinError(int(i[0]+1),int(i[3]+1),i[7])
h.GetXaxis().SetTitle("eta")
h.GetYaxis().SetTitle("pT")
h.Draw("COLZ")
h.SetOption("COLZ")

myFile = ROOT.TFile.Open("%s.root" %sys.argv[1], "RECREATE")
myFile.WriteObject(h,"%s" %sys.argv[2])
