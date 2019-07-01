#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:10:23 2018

@author: jmw254
"""

#import packages
import graph_tool.all as gt

#------------------------------------------------------
#Variables:
NETWORK_FILE=""
NETWORK_FEATURE_FILE=""
#------------------------------------------------------

#import the graph
g=gt.load_graph(NETWORK_FILE)
#calculate the features using inbuilt graph_tool functions

#Pagerank
rank=gt.pagerank(g)
print("pagerank has been calculated")

#HITS y-hubs and x-authorities
eigenvalue, xauthorities, yhubs =gt.hits(g)
print("HITS values have been calculated")

#betweenness centrality
between_vp,between_ep =gt.betweenness(g)
print("betweenness centrality has been calculated")

#save external to internal property map
#this makes the features accessible in the future when loading the graph
g.vertex_properties["page_rank"]=rank
g.vertex_properties["x_authorities"]=xauthorities
g.vertex_properties["y_hubs"]=yhubs
g.vertex_properties["b_cent"]=between_vp


g.save(NETWORK_FEATURE_FILE,fmt='gt')


