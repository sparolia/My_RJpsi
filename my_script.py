#macro that makes a stack plot with CMS requirements
import ROOT
import sys
from cmsstyle import CMS_lumi
#from sample import sample_dic
import os
import math
#import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
#from root_pandas import to_root
#from root_pandas import read_root
import random
var =["pt_miss_scal","m_miss_sq","Q_sq","pt_var"]
var_lowrange = [0,0,0,-30]
var_highrange = [30,15,15,30]
inFileName = sys.argv[1]
inFile = ROOT.TFile.Open(inFileName,"READ")
tree = inFile.Get("BTo3Mu")
#his = ROOT.TH1F("sample.histoName","MC" ,var.nbins,var.xmin,var.xmax)


for i in range(len(var)):
    his = ROOT.TH1F(var[i],"",100,var_lowrange[i],var_highrange[i])
    for entryNum in range (0,tree.GetEntries()):
        tree.GetEntry(entryNum)
        his.Fill(getattr (tree,var[i]))
        #his.Fill(pt_miss_scal)


    ROOT.gROOT.SetBatch()
    #ROOT.gStyle.SetOptStat(0)
    his.SetMarkerStyle(20)
    his.SetMarkerSize(0.9)
    his.SetLineColor(ROOT.kBlack)

    
    c = ROOT.TCanvas("","",700, 700)
    
    c.cd() 
    c.SetTicks(True)
    c.SetBottomMargin(2)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)
    c.Draw()
    c.cd()
    
    his.Draw("")
    #legend.Draw("SAME")
    #his.SetTitle("m^{2}_{miss}")
    his.GetYaxis().SetTitle("Events")
    his.GetYaxis().SetLabelOffset(0.01)
    his.GetYaxis().SetTitleOffset(2)
    his.GetYaxis().SetLabelFont(43)
    his.GetYaxis().SetLabelSize(15)
    his.GetYaxis().SetTitleFont(43)
    his.GetYaxis().SetTitleSize(18)
    his.SetMaximum(his.GetMaximum()*1.5)
    
    his.GetXaxis().SetLabelOffset(0.01)
    his.GetXaxis().SetTitleOffset(1.6)
    his.GetXaxis().SetLabelFont(43)
    his.GetXaxis().SetLabelSize(15)
    his.GetXaxis().SetTitleFont(43)
    his.GetXaxis().SetTitleSize(18)
    his.GetXaxis().SetTitle(var[i]+"(GeV)")
    CMS_lumi(c, 4, 1)
    c.RedrawAxis()
    c.SaveAs(var[i]+".png")
        
inFile.Close()
