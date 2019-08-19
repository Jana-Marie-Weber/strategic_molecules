#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:03:31 2019

@author: jmw254
"""

#NREL screening
#you nee dto specify the network you want to load
NETWORK_FILE=""


import cPickle as pickle
import graph_tool.all as gt

hubs=pickle.load(open("hubs_contamination0.001.p","rb"))
g=gt.load_graph(NETWORK_FILE)

NREL_BB=[]
with open("NREL_Building_Blocks.dat","rb") as infile:
    for line in infile:
        line=line.split(",")
        NREL_BB.append(line[0])
#----------------------------------------------------------------------
# Screening of NREL_BB in the network
#----------------------------------------------------------------------

NREL_BB_4HAP=[]
for v in g.vertices():
    if g.vp.xrn[v] in NREL_BB:
        NREL_BB_4HAP.append(g.vp.xrn[v])

#----------------------------------------------------------------------
# Screening of NREL_BB in the hubs
#----------------------------------------------------------------------

inter=set(hubs).intersection(set(NREL_BB))


#----------------------------------------------------------------------
# Saving results
#----------------------------------------------------------------------

with open("Results_NREL.txt","w") as outfile:
    outfile.write("Results from NREL analysis")
    outfile.write("\n")
    outfile.write("\n")
    outfile.write("We find: "+str(len(NREL_BB))+" in Reaxys.")
    outfile.write("\n")
    outfile.write("Out of these, we find: "+str(len(NREL_BB_4HAP))+" to be in the 4HAP network.")
    outfile.write("\n")
    outfile.write("Out of these: "+str(len(NREL_BB_4HAP))+" we find: "+str(len(inter))+ "classified as hubs by the algorithm.")
    outfile.write("\n")
    outfile.write("\n")
    outfile.write(str(inter))


