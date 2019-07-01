#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:11:07 2018

@author: jmw254
"""
#import all modules
from __future__ import division
import graph_tool.all as gt
import numpy as np
import cPickle as pickle




#load the limonene graph for screening of feedstock, hubs, and pharmaceuticals
g=gt.load_graph("limonene_cent_k_core.gt")
#get rid of parallel edges
gt.remove_parallel_edges(g)


#opens the smalles contamination rate tested. This will only evalutate the highest hubs found.
hubs_XRN=pickle.load(open("hubs_contamination0.0001.p","rb"))

#setting up lists for the hubs int value and rxn id
hubs=[]
hubs_xrn=[]

#setting status to 1 when performing the script the first time. Afterwards setting it to 2, so that we do not need to iterate through the whole network again.
status=2
if status==1:
    #list of all REAXIS IDs of top200 pharmaceuticals
    Top_200=[]
    
    # The output1.dat contains pharmaceutical products and their xrn id.
    with open("output1.dat") as file:
        #we need to split the names from the xrn id so that we can perform a search with a list of xrn ids.
        for line in file:
            splitline=line.split()
            RXID=splitline[1]
            Top_200.append(RXID)
            
    #we need to fill in int values of the hubs and the pharmaceuticals. Also, not all pharmaceuticals will be in the test network around limonene.
    in_network=[]
    in_network_xrn=[]        
    for v in g.vertices():
        if g.vp.xrn[v] in Top_200:
            in_network.append(int(v))
            in_network_xrn.append(g.vp.xrn[v])
        elif g.vp.xrn[v] in hubs_XRN:
            hubs.append(int(v))
            hubs_xrn.append(g.vp.xrn[v])

    #we save all of these lists for future use.
    pickle.dump(in_network,open("pharma_in_network.p","wb"))
    pickle.dump(in_network_xrn,open("pharma_in_network_xrn.p","wb"))
    pickle.dump(hubs,open("hubs_Int.p","wb"))
    pickle.dump(hubs_xrn,open("hubs_xrn.p","wb"))

#if the script isnt used for the first time and the status has been changed, than you will start here with loading the previously dumped files.
else:
    in_network=pickle.load(open("pharma_in_network.p","rb"))
    in_network_xrn=pickle.load(open("pharma_in_network_xrn.p","rb"))
    hubs=pickle.load(open("hubs_Int.p","rb"))
    hubs_xrn=pickle.load(open("hubs_xrn.p","rb"))



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#hardcoded biofeed input!These are the three Int values of the CST components
Biofeed=[9896190,345687,2239190]
#------------------------------------------------------------------------------
#you can specifiy in the bracets if you want to search over all hubs and for all pharmaceuticals, or if it should just be a fraction used for the search.
#in_network=in_network[]
#hubs=hubs[]

#the following are empty lists that will be filled when performing the DF search in the network.
paths_to_hub=[]
paths_from_hub=[]
none_zero_paths_from_hub=[]
variance=[]
Feedstocks=[]
Feedstocks_r=[]
#iteration through all hubs
for k in range(0,len(hubs)):
    #temp will store a list for each hubs and its connections to the biofeed.
    temp=[]
    #iteration through the feedstock components
    for j in range(0,len(Biofeed)):
        #a is set to 0 to start counting the number of possible paths from one feedstock to the hub.
        a=0
        #the next line is the DF search implemented from graph_tool. You need to specify here if you wish to change the screening range. we have used a cut-off value of 2 as we found the hubs closer to the feedstocks than to the pharmaceuticals.
        for path in gt.all_paths(g,Biofeed[j],hubs[k],cutoff=2):
            a=a+1
        print("biofeed:",j,"has:",a,"paths to the hub:",hubs_xrn[k])
        #the number off all paths from one feedtsocj to the hub is added to the temporary list.
        temp.append(a)
    #the temporary list is added to the overall storing list.
    paths_to_hub.append(temp)

    #in this piece of code we find how many feedstock components are connected in the screened range. For all feedstock components, we also check the maximum connectiviy between the hub and each single one.
    a_1=0
    limit=0
    for i in temp:
        if i!=0:
            a_1+=1
            maxi=max(i,limit)
            limit=maxi
    Feedstocks.append(a_1)
    Feedstocks_r.append(limit)
        
    
    #now we check the connections to the pharmaceuticals from the hub.
    paths=[]
    #iteration over all pharmaceuticals
    for i in range(0,len(in_network)):
        b=0
        #this is the piece of code that performs the DF search. If you wish to change the screening range, you need to do this here. It is highly recommended not to exceet the range 3, as the search drastically will slow down with any value greater than 3. 
        for path in gt.all_paths(g,hubs[k],in_network[i],cutoff=3):
            b=b+1
        #we add the findings to another temporary list (paths) and add this list to the overall storing list.
        paths.append(b)
    paths_from_hub.append(paths)

    #this piece of code checks how many of the pharmaceuticals are connected to the hub.
    c=0
    for i in paths:
        if i != 0:
            c=c+1            
    print("The hub",hubs_xrn[k], "has paths to:",c,"out of:",len(in_network),"pharmaceuticals")
    none_zero_paths_from_hub.append(c)

    #here we find the pair with the largest number of alternative connections.
    limit=0
    for i in paths:
        high=max(i,limit)
        limit=high
    variance.append(limit)

         

#relative numbers:
#we generate some relative numbers. The first and third ones are relative with reagrd to the overall number of feedstocks and pharmaceuticals respectively. The second and the forth element are relative to the maximal number of paths between pairs of all hubs. Hence, this relative value already gives an comparison between all hubs.
#1.
Feedstocks_relative=[x/len(Biofeed) for x in Feedstocks]
#2.
Feedstocks_r_relative=[x/max(Feedstocks_r) for x in Feedstocks_r]
#3.
Pharma_relative=[x/len(in_network) for x in none_zero_paths_from_hub]
#4.
Pharma_r_relative=[x/max(variance) for x in variance]

#print the results
print(Feedstocks_relative)
print(Feedstocks_r_relative)
print(Pharma_relative)
print(Pharma_r_relative)

#save the results for future use
pickle.dump(Feedstocks_relative,open("feeds.p","wb"))
pickle.dump(Feedstocks_r_relative,open("feeds_roues.p","wb"))
pickle.dump(Pharma_relative,open("pharma.p","wb"))
pickle.dump(Pharma_r_relative,open("pharma_roues.p","wb"))

#save the results in csv for plotting
evaluation_data=np.asarray([Feedstocks_relative,Feedstocks_r_relative,Pharma_relative,Pharma_r_relative])
np.savetxt("evaluation.csv", evaluation_data, delimiter=",")
np.savetxt("hubs_xrn.csv", hubs_xrn, delimiter=",", fmt='%s')
        
