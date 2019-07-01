#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 10:02:33 2019

@author: jmw254
"""

import __future__ 
#load all packages that you need
import graph_tool.all as gt
import sklearn
from sklearn import preprocessing
import numpy as np 
from sklearn.ensemble import IsolationForest
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import cPickle as pickle
import pandas as pd
import seaborn as sb
import matplotlib.cm as cm
import matplotlib.colors as colors

#------------------------------------------------------
#Input variables
NETWORK_FEATURE_FILE=""
#------------------------------------------------------


print(__doc__)


#load the graph and remove parallel links
g=gt.load_graph(NETWORK_FEATURE_FILE)
print("...removing parallel edges...")
gt.remove_parallel_edges(g)


#set up feature vectors
print("...setting up feature vector...")
#first we build empty lists for all features
a=[]
a_xrns=[]
cents=[]
pages=[]
hits=[]
degrees_o=[]
degrees_i=[]
xrn=[]
#second we iterate through the network filling the properties into the lists
for v in g.vertices():
    cents.append(g.vp.b_cent[v])
    pages.append(g.vp.page_rank[v])
    hits.append(g.vp.y_hubs[v])
    degrees_o.append(v.out_degree())
    degrees_i.append(v.in_degree())
    xrn.append(g.vp.xrn[v])
#third we describe one molecule by all features belong to it. We set up a matrix like list of lists in python.
for i in range(0,len(cents)):
    a.append([degrees_o[i],degrees_i[i],cents[i],pages[i],hits[i]])

    
#normalise the feature vector by max min normalisation from sklearn
print("...normalisation of feature vector...")    
min_max_scaler = sklearn.preprocessing.MinMaxScaler()
a_norm=min_max_scaler.fit_transform(a) 
   
#check the normalisation-------------------------------------------------------  
#------------------------------------------------------------------------------
#test the normalised feature range and print it to an external file for future references
maxis=np.amax(a_norm, axis=0)
minis=np.amin(a_norm, axis=0)

with open("normalisation_feature_vector.txt","w") as outfile:
    outfile.write("outdegrees max: " + str(max(degrees_o))+ ", min: " +str(min(degrees_o)))
    outfile.write("\n")
    outfile.write("indegrees max :" + str(max(degrees_i)) +", min: " + str(min(degrees_i)))
    outfile.write("\n")
    outfile.write("cents max: "+ str(max(cents))+ ", min: "+ str(min(cents)))
    outfile.write("\n")
    outfile.write("pages max: "+ str(max(pages)) +", min: "+ str(min(pages)))
    outfile.write("\n")
    outfile.write("hits max: "+ str(max(hits))+ ", min: "+ str(min(hits)))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("new out degree max: "+ str(maxis[0])+ ", new out degree min: "+ str(minis[0]))
    outfile.write("\n")
    outfile.write("new in degree max: "+ str(maxis[1])+ ", new in degree min: "+ str(minis[1]))
    outfile.write("\n")
    outfile.write("new betweenness max: "+str(maxis[2])+ ", new betweenness min: "+ str(minis[2]))
    outfile.write("\n")
    outfile.write("new pages max: "+ str(maxis[3])+ ", new pages min: "+ str(minis[3]))
    outfile.write("\n")
    outfile.write("new hits max: " + str(maxis[4]) + ", new hits min: "+ str(minis[4]))
    outfile.write("\n")
#------------------------------------------------------------------------------------
    
#----------------------------------------------------------------------
#perform PCA on the data.
#we perform PCA two times. The first one shows all principal components and their relevance. The second one then reduces to two principal components.
#we write the results to a text file for future references.
print("...performing PCA...")
pca = PCA(n_components=5)
X_all_compoenents_new=pca.fit_transform(a_norm)
PCA(copy=True, iterated_power='auto', n_components=5, random_state=None,
    svd_solver='auto', tol=0.0, whiten=False)

#write PCA results to file
with open("PCA_results_all_components.txt","w") as outfile:
    outfile.write("this is the explained_variance: ")
    outfile.write(str(pca.explained_variance_))
    outfile.write('\n')
    outfile.write("this is the explained_variance ratio: ")
    outfile.write(str(pca.explained_variance_ratio_))
    outfile.write('\n')
    outfile.write("this is the singular value: ")
    outfile.write(str(pca.singular_values_))
    outfile.write("\n")
    outfile.write("this is the explanined variance sum:")
    outfile.write(str(pca.explained_variance_ratio_.sum()))

#we wish to investigating how much each feature contributes to which principal component.
#Hence, in the order of the a matrix, which contain all the data, we name the columns.
variable_names=['out degree', 'in degree', 'betweenness', 'pagerank', 'hits y-hubs']
comps=pd.DataFrame(pca.components_,columns=variable_names)
print(comps)

#generate a heat map of the feature contribution to the principal components
add=1
if add==1:
    sb.heatmap(comps)
    plt.savefig("heatmap_correlation_all.pdf") 
    
    
#second performance of PCA to compress the data. 
pca = PCA(n_components=2)
Xnew=pca.fit_transform(a_norm)
PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,
    svd_solver='auto', tol=0.0, whiten=False)

#write PCA results to file
with open("PCA_results_2_components.txt","w") as outfile:
    outfile.write("this is the explained_variance: ")
    outfile.write(str(pca.explained_variance_))
    outfile.write('\n')
    outfile.write("this is the explained_variance ratio: ")
    outfile.write(str(pca.explained_variance_ratio_))
    outfile.write('\n')
    outfile.write("this is the singular value: ")
    outfile.write(str(pca.singular_values_)) 
    outfile.write("\n")
    outfile.write("this is the explanined variance sum:")
    outfile.write(str(pca.explained_variance_ratio_.sum()))
    
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#outlier detection with the isolation forest  

# these contamination rates will be tested
contamination_rate=[0.002,0.0015,0.001,0.00095,0.0009,0.00085,0.0008,0.00075,0.0007,0.00065,0.0006,0.00055,0.0005]
#contamination_rate=[0.0001]
#this dictionary will hold teh number of hubs classified by each of the contamination rates
hubs_per_rate={}

#the routine begins here:
for rate in contamination_rate:
    # fit the model, all samples are used for all trees. No sampling takes place.
    clf = IsolationForest(max_samples=len(Xnew), contamination=rate)
    clf.fit(Xnew)
    #a prediction for all samples is performed. The value for outliers will be -1.
    y_pred = clf.predict(Xnew) 
    
    #transformation of the PCA-based feature array into a numpy array for multiple sampling
    xnew=np.array(Xnew)
    
    #lists for hubs and non hubs and their PC1 and PC2 are set up. hubs_1 is the value of PC1 for all moelcules classified as hubs. hubs_2 is the value of PC2 for all molecules classified as hubs. The general hubs list will be filled up with the xrns of all hubs.
    hubs_1=[]
    hubs_2=[]
    non_hubs_1=[]
    non_hubs_2=[]
    hubs=[]
    #An iteration through the prediction vector allows to sort the moelcules.
    for i in range(0,len(y_pred)):
        if y_pred[i]==-1:
            hubs.append(xrn[i])
            hubs_1.append(xnew[(i,0)])
            hubs_2.append(xnew[(i,1)])
        else:
            non_hubs_1.append(xnew[(i,0)])
            non_hubs_2.append(xnew[(i,1)])
    #these lists are now saved to csv format to use them in origin
    np.savetxt("hubs_PC0_rate"+str(rate)+".csv", hubs_1, delimiter=",", fmt='%s')
    np.savetxt("hubs_PC1_rate"+str(rate)+".csv", hubs_2, delimiter=",", fmt='%s')
    np.savetxt("non_hubs_PC0_rate"+str(rate)+".csv", non_hubs_1, delimiter=",", fmt='%s')
    np.savetxt("non_hubs_PC1_rate"+str(rate)+".csv", non_hubs_2, delimiter=",", fmt='%s')
    np.savetxt("hubs"+str(rate)+".csv", hubs, delimiter=",", fmt='%s')



    #the results are then plotted in jpg for a rough evaluation
    plt.clf()
    plt.plot(hubs_1,hubs_2,'bo',non_hubs_1,non_hubs_2,'ro', markersize=3)
    plt.axis([-0.001,0.03,-0.02,0.005])
    #plt.axis([min(PC1)-0.2,max(PC1)+0.2,min(PC2)-0.2,max(PC2)+0.2])
    plt.xlabel("Principal component 0",fontsize=22)
    plt.ylabel("Principal component 1",fontsize=22)
    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(bottom=0.2)
    plt.tick_params(axis="both", which='major', labelsize=14)
    plt.legend(("classified as outliers","classified as normal"),
        loc='lower left', frameon=False, fontsize= 18)
    plt.savefig("contamination_rate"+str(rate)+".jpg")

    #the hubs xrns are saved for future use with pickle
    pickle.dump(hubs,open("hubs_contamination"+str(rate)+".p", "wb"))
    #the number of all hubs with this contamination rate is displayed.
    print(len(hubs))
    #this number is added to the overall dictionary
    hubs_per_rate[rate]=len(hubs)
    
    
with open("number_hubs_per_rate.txt","w") as outfile:
    outfile.write(str(hubs_per_rate))








    
