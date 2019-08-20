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

#-------------------------------------
#Variables
SCREENING_NETWORK=""
PHARMA_LIST="pharmaceuticals.dat"
#------------------------------------

#loading of a second graph used in previous studies. It is larger and shows more than 90% of the 4HAP graph. It potentially includes more feedstock and pharmaceutical regions.
g=gt.load_graph(SCREENING_NETWORK)
print("...removing parallel edges...")
gt.remove_parallel_edges(g)


#loading of the pharmalist and the hubs---------------------------------------
#either loaded the first time or just loaded from pickle----------------------

#if you run the same script again, you can use the pickled versions of the lists instead of iterating through everything again.
first_time=1
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
    with open(PHARMA_LIST) as file:
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
#-----------Evaluation from CST to hubs----------------------------------------
#These are the Int Values of the three CST components in the limonene network! BE CAREFUL! when using another network input the xrns and find the Int values!!
Biofeed=[9896190,345687,2239190]

#evaluation steps with cutoff 1 are performed
print("...evaluation ways to hubs...cutoff= 1")

#setting up an empty list for all hubs that are reached from the biofeedstock components
hubs_reached=[]
#iteration through all hubs
for k in range(0,len(hubs)):
    #in this list, we save the number of path between the feedstcok components and the hub. It will look like e.g. temp=[1,0,6]
    temp=[]
    #iteration through all biofeeds components
    for j in range(0,len(Biofeed)):
        #for each component we want to start counting the paths from new
        a=0
        #path search inbuilt in graph_tools. Starting point is the biofeedstock and end point is the hub. Only one reaction step is allowed.
        for path in gt.all_paths(g,Biofeed[j],hubs[k],cutoff=1):
            #all possible paths are counted
            a=a+1
        #the temp list is filled with the number of paths for one hub and all biofeeds
        temp.append(a)
    #the temp list is added to an overall list, so that it can be emptied for the next hub.
    hubs_reached.append(temp)

print(hubs_reached)

#an dictionary is set to count all hubs that are conntected to the feedstock.
connection_to_hubs={}
#iteration through all hubs
for i in range(0,len(hubs_reached)):
    for k in hubs_reached[i]:
        #if any of the paths from the feedstocks to one hub counts more than 0, this molecule will be regarded for the dictionary. If more than one of the paths is larger than 0 the dictionary entry will be overwritten without any loss as it is the same key and the same value.
        if not k ==0:
            connection_to_hubs[hubs_xrn[i]]=hubs_reached[i]

            
            

#a negative probe is set up using randomly drawn molecules.
negative_probe=[]
negative_probe_xrn=[]
#a subset of the same length as the hubs is assembled
while len(negative_probe) < len(hubs):
    temp=random.randint(0,len(all_ids))
    #we prevent duplicate draws out of the overall set
    if temp not in negative_probe:
        negative_probe.append(all_ids[temp])
        negative_probe_xrn.append(all_xrns[temp])

#now the same evaluation as to the hubs will take place to the random subsets
print("...evaluating ways to negative_probe...cutoff= 1")
       
#see all comments above. This is the same piece of code but for the negative probe instead of the hubs.
negative_hubs_reached=[]
for k in range(0,len(negative_probe)):
    temp=[]
    for j in range(0,len(Biofeed)):
        a=0
        for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=1):
            a=a+1
        temp.append(a)
    negative_hubs_reached.append(temp)

print(negative_hubs_reached)

connection_to_negative={}
for i in range(0,len(negative_hubs_reached)):
    for k in negative_hubs_reached[i]:
        if not k ==0:
            connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]

#both results will be recorded in a text file. We record for hubs and the random sample: The cut-off value, the overall number of connected molecules, the dictionary values for all three feedstcok components.
with open("evaluation1.txt","w") as outfile:
    outfile.write("CST to hubs in cutoff 1:-----------------") 
    outfile.write("\n") 
    outfile.write("number of hubs reached: "+ str(len(connection_to_hubs)))
    outfile.write("\n") 
    outfile.write(str(connection_to_hubs))
    outfile.write("\n")
    outfile.write("CST to random in cutoff 1:--------------------")
    outfile.write("\n")
    outfile.write("number of randoms reached: "+str(len(connection_to_negative))) 
    outfile.write("\n")
    outfile.write(str(connection_to_negative))
    
 #-----------------------------------------------------------------------------   

    
#we perfrom excatly the same step as before for a cut-off value of 2. Please see the documentation of the above code.
print("...evaluation ways to hubs...cutoff= 2")    



hubs_reached=[]
for k in range(0,len(hubs)):
    temp=[]
    for j in range(0,len(Biofeed)):
        a=0
        for path in gt.all_paths(g,Biofeed[j],hubs[k],cutoff=2):
            a=a+1
        temp.append(a)
    hubs_reached.append(temp)

print(hubs_reached)

connection_to_hubs={}
for i in range(0,len(hubs_reached)):
    for k in hubs_reached[i]:
        if not k ==0:
            connection_to_hubs[hubs_xrn[i]]=hubs_reached[i]

            
            


negative_probe=[]
negative_probe_xrn=[]
while len(negative_probe) < len(hubs):
    temp=random.randint(0,len(all_ids))
    if temp not in negative_probe:
        negative_probe.append(all_ids[temp])
        negative_probe_xrn.append(all_xrns[temp])
 
print("...evaluating ways to negative_probe...cutoff= 2")
       
        
negative_hubs_reached=[]
for k in range(0,len(negative_probe)):
    temp=[]
    for j in range(0,len(Biofeed)):
        a=0
        for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=2):
            a=a+1
        temp.append(a)
    negative_hubs_reached.append(temp)

print(negative_hubs_reached)

connection_to_negative={}
for i in range(0,len(negative_hubs_reached)):
    for k in negative_hubs_reached[i]:
        if not k ==0:
            connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]


with open("evaluation2.txt","w") as outfile:
    outfile.write("CST to hubs in cutoff 2:-----------------") 
    outfile.write("\n") 
    outfile.write("number of hubs reached: "+ str(len(connection_to_hubs)))
    outfile.write("\n") 
    outfile.write(str(connection_to_hubs))
    outfile.write("\n")
    outfile.write("CST to random in cutoff 2:--------------------")
    outfile.write("\n")
    outfile.write("number of randoms reached: "+str(len(connection_to_negative))) 
    outfile.write("\n")
    outfile.write(str(connection_to_negative))
    

 #-----------------------------------------------------------------------------

#we perfrom excatly the same step as before for a cut-off value of 3. Please see the documentation of the above code.
print("...evaluation ways to hubs...cutoff= 3")    


hubs_reached=[]
for k in range(0,len(hubs)):
    temp=[]
    for j in range(0,len(Biofeed)):
        a=0
        for path in gt.all_paths(g,Biofeed[j],hubs[k],cutoff=3):
            a=a+1
        temp.append(a)
    hubs_reached.append(temp)

print(hubs_reached)

connection_to_hubs={}
for i in range(0,len(hubs_reached)):
    for k in hubs_reached[i]:
        if not k ==0:
            connection_to_hubs[hubs_xrn[i]]=hubs_reached[i]

            
            


negative_probe=[]
negative_probe_xrn=[]
while len(negative_probe) < len(hubs):
    temp=random.randint(0,len(all_ids))
    if temp not in negative_probe:
        negative_probe.append(all_ids[temp])
        negative_probe_xrn.append(all_xrns[temp])
 
print("...evaluating ways to negative_probe...cutoff= 3")
       
        
negative_hubs_reached=[]
for k in range(0,len(negative_probe)):
    temp=[]
    for j in range(0,len(Biofeed)):
        a=0
        for path in gt.all_paths(g,Biofeed[j],negative_probe[k],cutoff=3):
            a=a+1
        temp.append(a)
    negative_hubs_reached.append(temp)

print(negative_hubs_reached)

connection_to_negative={}
for i in range(0,len(negative_hubs_reached)):
    for k in negative_hubs_reached[i]:
        if not k ==0:
            connection_to_negative[negative_probe_xrn[i]]=negative_hubs_reached[i]


with open("evaluation3.txt","w") as outfile:
    outfile.write("CST to hubs in cutoff 3:-----------------") 
    outfile.write("\n") 
    outfile.write("number of hubs reached: "+ str(len(connection_to_hubs)))
    outfile.write("\n") 
    outfile.write(str(connection_to_hubs))
    outfile.write("\n")
    outfile.write("CST to random in cutoff 3:--------------------")
    outfile.write("\n")
    outfile.write("number of randoms reached: "+str(len(connection_to_negative))) 
    outfile.write("\n")
    outfile.write(str(connection_to_negative))    



                


