# Identification of strategic molecules

This repository contains the code for the network-based pipeline for the identification of strategic molecules by [(Weber et al., 2019)]. Due to license agreements the datamining part is not part of this repository. The pipeline outlines steps 3 to 5 assuming a network-based data collection had been assembled previously. To cite from this work use [(Weber et al., 2019)]. 




<img align="centre" src="documents/pipeline.png" width="600" > 

## Requirements
We used python 2.7.12 (updated Nov 2018) for the computation. System configurations can be seen in the [requirment.txt] file. Graph-tool was installed from source, as it has C++ dependencies, for future use we would advise to use the graph-tool docker as recommended on the [graph-tool webpage]. 

## Using the algorithm 

To use the strategic molecules in **a similar region of chemistry**, simply refer to [(Weber et al., 2019)] and use the full list of strategic molecules including structures and Reaxys IDs from the Electronic Supplementary Information (ESI). This data set is centred around 4-hydroxyacetophenone and contains up to three reaction steps from the centre.

To use this algorithm on your own set of **chemical reaction data**, build a network structure based on your data, e.g. using [graph-tool]. Molecules are represented by nodes and reactions are represented by edges. We wired all reactants to all products and only attached molecules identification numbers to the nodes using graph-tools property maps [(Jacob et al., 2018)]. After network assembly, simply download all files contained in this repository and start the pipeline with the feature computation. The variable `NETWORK_FILE` in [feature_calculation.py] should be adjusted to your network in ".gt" format and the variable `NETWORK_FEATURE_FILE` to the name you wish to save the network with features as. 

The general workflow of the pipeline can be applied to **different network types** as well. Build a network structure based on your data, e.g. using [graph-tool], and download all files from this repository. You might need to manually modify the selected features so that they best describe optimal locations in your network problem. Before applying the isolation forest outlier detection part, verify that optimal locations rank exceptional in your feature value distributions. 

## Repository organisation 
The folders 1_feature_computation, 2_isolation_forest, and 3_screening each contain parts of the pipeline and their own README.md explaining code and implementation. Please specify the variable `NETWORK_FILE` in [feature_calculation.py] in 3_feature_computation and follow the instructions throughout the rest of the pipeline. For questions please contact: <jmw254@cam.ac.uk>.



[(Weber et al., 2019)]: https://pubs.rsc.org/en/content/articlelanding/2019/re/c9re00213h/unauth#!divAbstract
[graph-tool]: https://graph-tool.skewed.de/static/doc/quickstart.html
[feature_calculation.py]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/1_feature_computation/feature_calculation.py
[(Jacob et al., 2018)]: https://pubs.rsc.org/en/content/articlehtml/2018/re/c7re00129k

[graph-tool webpage]: https://git.skewed.de/count0/graph-tool/wikis/installation-instructions#installing-using-docker
[requirment.txt]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/requirements.txt
