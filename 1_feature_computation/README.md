# Graph features

In this part of the work, we need to compute the features to describe important molecules (optimal locations) in networks. For the purpose of this work we considered both [in and out degree of nodes], [betweenness centrality], [HITS y-hubs], and the [pagerank] of a node. These were chosen based on a literature review documented in [(Weber et al., 2019)]. 

For reproducibility of other network problems, an intensive literature review should be performed before chosing the descriptive features. Here, for important molecules we consider following aspects: 

- participation in chemical reactions (shown by in and out degree of a node),

- long-range importance, e.g. importance of molecule in multi-step reactions (betweenness centrality), and

- closeness to often used molecules/ importnat industrial intermediates (pagerank and HITS y-hubs). 

In the [feature_calculation.py] script we use the [graph-tools] library to compute all features. Starting point is a graph-based data structure saved in the ".gt" format. The only information attached to each node at this state of the pipeline is the Reaxys ID. The values are then saved to each molecules in the format of a property map, making them accessible each time one opens the newly saved graph. 




[in and out degree of nodes]:https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.Vertex.out_degree
[betweenness centrality]: https://graph-tool.skewed.de/static/doc/centrality.html?highlight=betweenness#graph_tool.centrality.betweenness
[HITS y-hubs]: https://graph-tool.skewed.de/static/doc/centrality.html?highlight=betweenness#graph_tool.centrality.hits
[pagerank]: https://graph-tool.skewed.de/static/doc/centrality.html?highlight=betweenness#graph_tool.centrality.pagerank
[(Weber et al., 2019)]: https://pubs.rsc.org/en/content/articlehtml/2018/re/c7re00129k
[graph-tools]: https://graph-tool.skewed.de
[feature_calculation.py]:https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/1_feature_computation/feature_calculation.py
