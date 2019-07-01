#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 10:43:57 2019

@author: jmw254
"""
#load all modules that you need
import __future__
import graph_tool.all as gt
import numpy as np 
import cPickle as pickle
import random 


#loading of a second graph used in previous studies. It is larger and shows more than 90% of the 4HAP graph. It potentially includes more feedstock and pharmaceutical regions.
g=gt.load_graph("limonene_cent_k_core.gt")
print("...removing parallel edges...")
gt.remove_parallel_edges(g)


#loading of the pharmalist and the hubs---------------------------------------
#either loaded the first time or just loaded from pickle----------------------

#if you run the same script again, you can use the pickled versions of the lists instead of iterating through everything again.
first_time=2
if first_time==1:
    #setting up empty lists for the pharmaceuticals. We need both, one list saveing the xrn ID and another list saving the Int value of the node in the network.
    Top_200=[]   
    T200_in_network=[]
    T200_in_network_xrn=[]
    #setting up empty lists for the hubs. We need both, one list saveing the xrn ID and another list saving the Int value of the node in the network.
    hubs_XRN=pickle.load(open( "hubs_contamination0.001.p", "rb" ))
    hubs=[]
    hubs_xrn=[]
    #setting up empty lists for the whole network. These will be used to draw random samples from the network
    all_ids=[]
    all_xrns=[]
    
    #opening of the pharmaceutical inputs and saving their xrn ID.
    with open("output1.dat") as file:
        for line in file:
            splitline=line.split()
            RXID=splitline[1]
            Top_200.append(RXID)

    #iterating through the whole network and finding both the pharmaceuticals and hubs. The int value of these in the network needs to be saved.
    for v in g.vertices():
        if g.vp.xrn[v] in Top_200:
            T200_in_network.append(int(v))
            T200_in_network_xrn.append(g.vp.xrn[v])
        elif g.vp.xrn[v] in hubs_XRN:
            hubs.append(int(v))
            hubs_xrn.append(g.vp.xrn[v])
        #the same is done for every single molecule in the network. Here we will draw the random subsets from. 
        all_ids.append(int(v))
        all_xrns.append(g.vp.xrn[v])
            
    #all assembled lists are pickled for future use
    pickle.dump(T200_in_network,open("pharma_in_network.p","wb"))
    pickle.dump(T200_in_network_xrn,open("pharma_in_network_xrn.p","wb"))
    pickle.dump(hubs,open("hubs_Int.p","wb"))
    pickle.dump(hubs_xrn,open("hubs_xrn.p","wb"))
    pickle.dump(all_ids,open("ids_of_all.p","wb"))
    pickle.dump(all_xrns,open("xrns_of_all.p","wb"))

#if running the script a second time, use the pickled version to save CP time as iterating through the network and comparing xrn Ids is a major operation.
else:

    T200_in_network=pickle.load(open("pharma_in_network.p","rb"))
    T200_in_network_xrn=pickle.load(open("pharma_in_network_xrn.p","rb"))
    hubs=pickle.load(open("hubs_Int.p","rb"))      
    hubs_xrn=pickle.load(open("hubs_xrn.p","rb")) 
    all_ids=pickle.load(open("ids_of_all.p","rb"))
    all_xrns=pickle.load(open("xrns_of_all.p","rb"))

print("...data is loaded...")


#------------------------------------------------------------------------------
#-----------Evaluation from CST to random subsets----------------------------------------
#These are the Int Values of the three CST components in the limonene network! BE CAREFUL! when using another network input the xrns and find the Int values!!
Biofeed=[9896190,345687,2239190]

#negative probes are set up using randomly drawn molecules. We set up 10 probes and use these for all cut-off values.
all_negative_probes=[]
all_negative_probes_xrns=[]
for h in range(0,10):
    negative_probe=[]
    negative_probe_xrn=[]
    #a subset of the same length as the hubs is assembled
    while len(negative_probe) < len(hubs):
        temp=random.randint(0,len(all_ids))
        #we prevent duplicate draws out of the overall set
        if temp not in negative_probe:
            negative_probe.append(all_ids[temp])
            negative_probe_xrn.append(all_xrns[temp])
    all_negative_probes.append(negative_probe)
    all_negative_probes_xrns.append(negative_probe_xrn)



#now the same evaluation as to the hubs will take place to the random subsets
print("...evaluating ways to negative_probe...cutoff= 1")
       
#see all comments for the hub molecule search. This is the same piece of code as used for the hubs before. Now we evaluate the random samples with a cut-off 1.
all_neg_hubs=[]
all_neg_dict=[]
for neg in range(0,len(all_negative_probes)):
    negative_hubs_reached=[]
    negative_probe=all_negative_probes[neg]
    for k in range(0,len(negative_probe)):
        temp=[]
        for j in range(0,len(Biofeed)):
            a=0
            for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=1):
                a=a+1
            temp.append(a)
        negative_hubs_reached.append(temp)
    all_neg_hubs.append(negative_hubs_reached)


    connection_to_negative={}
    for i in range(0,len(negative_hubs_reached)):
        for k in negative_hubs_reached[i]:
            if not k ==0:
                connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]
    all_neg_dict.append(connection_to_negative)


#The following will be recorded: The cut-off value, the overall number of connected molecules, the dictionary values for all three feedstcok components.
with open("evaluation_cut_off_1.txt","w") as outfile:
    outfile.write("CST to random in cutoff 1:--------------------")
    outfile.write("\n")
    outfile.write("------------------------------------------------------------------------------------")
    outfile.write("\n")
    for i in range(0,len(all_neg_hubs)):
        outfile.write("number of randoms reached in set " +str(i)+" is:"+str(len(all_neg_dict[i])))
        outfile.write("\n")
        outfile.write(str(all_neg_dict[i]))
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
print("...evaluating ways to negative_probe...cutoff= 2")

#see all comments above. This is the same piece of code for cut-off 2.
all_neg_hubs=[]
all_neg_dict=[]
for neg in range(0,len(all_negative_probes)):
    negative_hubs_reached=[]
    negative_probe=all_negative_probes[neg]
    for k in range(0,len(negative_probe)):
        temp=[]
        for j in range(0,len(Biofeed)):
            a=0
            for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=2):
                a=a+1
            temp.append(a)
        negative_hubs_reached.append(temp)
    all_neg_hubs.append(negative_hubs_reached)


    connection_to_negative={}
    for i in range(0,len(negative_hubs_reached)):
        for k in negative_hubs_reached[i]:
            if not k ==0:
                connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]
    all_neg_dict.append(connection_to_negative)


#The following will be recorded: The cut-off value, the overall number of connected molecules, the dictionary values for all three feedstcok components.
with open("evaluation_cut_off_2.txt","w") as outfile:
    outfile.write("CST to random in cutoff 2:--------------------")
    outfile.write("\n")
    outfile.write("------------------------------------------------------------------------------------")
    outfile.write("\n")
    for i in range(0,len(all_neg_hubs)):
        outfile.write("number of randoms reached in set " +str(i)+" is:"+str(len(all_neg_dict[i])))
        outfile.write("\n")
        outfile.write(str(all_neg_dict[i]))
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
print("...evaluating ways to negative_probe...cutoff= 3")

#see all comments above. This is the same piece of code but for cut-off 3.
all_neg_hubs=[]
all_neg_dict=[]
for neg in range(0,len(all_negative_probes)):
    negative_hubs_reached=[]
    negative_probe=all_negative_probes[neg]
    for k in range(0,len(negative_probe)):
        temp=[]
        for j in range(0,len(Biofeed)):
            a=0
            for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=3):
                a=a+1
            temp.append(a)
        negative_hubs_reached.append(temp)
    all_neg_hubs.append(negative_hubs_reached)


    connection_to_negative={}
    for i in range(0,len(negative_hubs_reached)):
        for k in negative_hubs_reached[i]:
            if not k ==0:
                connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]
    all_neg_dict.append(connection_to_negative)


#The following will be recorded: The cut-off value, the overall number of connected molecules, the dictionary values for all three feedstcok components.
with open("evaluation_cut_off_3.txt","w") as outfile:
    outfile.write("CST to random in cutoff 3:--------------------")
    outfile.write("\n")
    outfile.write("------------------------------------------------------------------------------------")
    outfile.write("\n")
    for i in range(0,len(all_neg_hubs)):
        outfile.write("number of randoms reached in set " +str(i)+" is:"+str(len(all_neg_dict[i])))
        outfile.write("\n")
        outfile.write(str(all_neg_dict[i]))
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")
        outfile.write("------------------------------------------------------------------------------------")
        outfile.write("\n")




