# Graph features

In this part of the work, we need to compute the features to describe important molecules (optimal locations) in networks. For the purpose of this work we considered both in and out degree of nodes, betweenness centrality, HITS y-hubs, and the pagerank of a node. These were chosen based on a literature review documented in (Weber et al., 2019). 

For reproducibility of other network problems, an intensive literature review should be performed before chosing the descriptive features. Here, for important molecules we consider following aspects: 

- participation in chemical reactions (shown by in and out degree of a node),

- long-range importance, e.g. importance of molecule in multi-step reactions (betweenness centrality), and

- closeness to often used molecules/ importnat industrial intermediates (pagerank and HITS y-hubs). 

In the feature_calculation.py script we use the inbuilt graph_tools library to compute all features. Starting point is a graph-based data structure saved in the .gt format. The only information attached to each node at this state of the pipeline is the Reaxys ID. The values are then saved to each molecules in the format of a property map, making them accessible each time one opens the newly saved graph. 

