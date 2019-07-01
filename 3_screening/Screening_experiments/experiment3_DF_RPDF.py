#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 09 11:45:34 2019

@author: jmw254
"""

from __future__ import division
import graph_tool.all as gt
import numpy as np
import cPickle as pickle
import time
import sys

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#INPUT
element=int(sys.argv[1])
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
print(element)

#the limonene graph is loaded into the system and parallel edges are removed.
g=gt.load_graph("limonene_cent_k_core.gt")
gt.remove_parallel_edges(g)

#a timer is started for the CPU comparison of the two methods
start=time.clock()


in_network=pickle.load(open("pharma_in_network.p","rb"))
in_network_xrn=pickle.load(open("pharma_in_network_xrn.p","rb"))
HUBS_ID=pickle.load(open("hubs_Int.p","rb"))
HUBS_XRN=pickle.load(open("hubs_xrn.p","rb"))

# double-check if the loaded files are correct
print("testing:There are:", len(HUBS_ID), len(HUBS_XRN)," hubs")



#hardcoded biofeed input! change as soon as mre data is available
Biofeed=[9896190,345687,2239190]


Feed=Biofeed[0]
Pharma=in_network[element]
Pharma_xrn=in_network_xrn[element]

#--------------------------------------------------------------------------
#This is the normal DF screening part
#--------------------------------------------------------------------------
easy_screen_3=[]
easy_screen_4=[]
easy_screen_5=[]
check10=time.clock()
for path in gt.all_paths(g,Feed,Pharma,cutoff=3):
    easy_screen_3.append(path)
print(len(easy_screen_3))
check11=time.clock()
print("time to screen with cutoff 3:", check11-check10)

check12=time.clock()
for path in gt.all_paths(g,Feed,Pharma,cutoff=4):
    easy_screen_4.append(path)
print(len(easy_screen_4))
check13=time.clock()
print("time to screen with cutoff 4:", check13-check12)
    


check14=time.clock()
for path in gt.all_paths(g,Feed,Pharma,cutoff=5):
    easy_screen_5.append(path)
print(len(easy_screen_5))
check15=time.clock()
print("time to screen with cutoff 5:", check15-check14)
#--------------------------------------------------------------------------
#This is the RPDF screening part
#--------------------------------------------------------------------------
paths_to_hub=[]
paths_to_hub2=[]
paths_from_hub=[]
paths_from_hub2=[]
none_zero_paths_from_hub=[]
variance=[]
Feedstocks=[]
Feedstocks_r=[]
hubs=HUBS_ID
    
check2=time.clock()
for k in range(0,len(hubs)):
    temp=[]
    temp_path=[]
    temp_p=[]
    for path in gt.all_paths(g,Biofeed[0],hubs[k],cutoff=2):
        temp_p.append(path)
    temp.append(len(temp_p))
    temp_path.append(temp_p)
    paths_to_hub.append(temp)
    paths_to_hub2.append(temp_path)
        
    paths=[]
    paths_one_hub_all_ends=[]
    temp_to_one_end_product=[]
    for path in gt.all_paths(g,hubs[k],Pharma,cutoff=3):
        temp_to_one_end_product.append(path)
    paths.append(len(temp_to_one_end_product))
    paths_one_hub_all_ends.append(temp_to_one_end_product)
    paths_from_hub.append(paths)
    paths_from_hub2.append(paths_one_hub_all_ends)

    
general=[]
for i in range(0,len(paths_to_hub)):
    temp1=paths_to_hub[i]
    temp2=paths_from_hub[i]
    for k in range(0,len(temp1)):
        temp3=temp1[k]*temp2[k]
        general.append(temp3)
all_combined_paths=sum(general)
check3=time.clock()





#--------------------------------------------------------------------------
#Saving the results
#--------------------------------------------------------------------------

CPU_Times=np.asarray([check11-check10,check13-check12,check15-check14])
np.savetxt("DF_cpu_times_"+str(Pharma_xrn)+".csv", CPU_Times, delimiter=",")

Paths=np.asarray([len(easy_screen_3),len(easy_screen_4),len(easy_screen_5)])
np.savetxt("DF_paths_"+str(Pharma_xrn)+".csv", Paths, delimiter=",")

CPU_Times_RPDF=np.asarray([check3-check2,all_combined_paths])
np.savetxt("RPDF_1_CPU_2_paths_"+str(Pharma_xrn)+".csv", CPU_Times_RPDF, delimiter=",")


