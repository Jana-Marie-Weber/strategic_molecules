# Screening 

Different screening methods are performed on basis of strategic moelcules. As outlined in [(Weber et al., 2019)] the follwoing experiments were conducted: 

1. Experiment no. 1: Depth-first (DF) searches from a renewable feedstock component to sets of randomly chosen molecules and to the set of strategic molecules are performed. 
2. Reduced Paths Depth First (RPDF) searches are conducted over strategic molecules to compare them with each other. Three steps from CST to strategic molecules and three steps from the strategic molecules to pharmaceuticals are tested.
3. DF and RPDF screening are compared with regard to the number of paths and CPU times needed for each screening.


## Experiment no. 1

The script [experiment1_strategic_molecules_random.py] tests all the cut-offs for the set of strategic moelcules and one set of random subsets. 

The script [experiment1_random.py] is used to test 10 sets of random moelcules within each cut-off.

The first part of the scripts loads the network used for screening `SCREENING_NETWORK`, the list of strategic molecules used `hubs_contamination0.001.p`, and a list of Top200 pharmaceuticals `output1.dat`. The pharmaceuticals are used in the following experiments. 

The second part of the script defines the three considered biologial feedstock components and assembles sets of randomly chosen molecules in the network. It is important to note that the DF screening algorithm works on the [Index number (numbering within graph-tools)] of the nodes and not on the given names. Hence, index numbering of nodes for the feedstock components were looked up previously and hardcoded here. 

The third part of the script performs the DF search of within cut-off value 1, 2 or 3 from the feedstock components to the sets of random subsets, `negative_probes`, and/or the strategic molecules and saves the results. 

## Experiment no. 2
The script [experiment2_strategic_molecuels_ranking.py] evaluates the connectivity of the strategic molecules within reaction paths from feedstock crude sulphate turpentine (CST) components to pharmaceutical products. 

The algorithm operates on the highest 56 strategic molecules, hence relative ranking values are computed between these. 

1. The network, the pharmaceutical components, and the strategic molecules are loaded.
2. DF search from feedstocks to strategic molecules is performed.
3. DF search from strategic molecules to pharmaceutical end products is perfomed.
4. Relative connectivity values are computed for all tested strategic molecules.
5. The results are saved. 


## Experiment no. 3
For this experiment, we perform DF and RPDF screening with each other. The starting point for all searches are the CST components, desired pharmaceutical products are inputed via the command line through sys input. 

The script [experiment3_DF_RPDF.py] computes CPU times and the number of path for the following searches:

- DF cut-off 3
- DF cut-off 4
- DF cut-off 5
- RPDF cut-off 5

The script [experiment3_RPDF.py] computes CPU times and the number rof path for RPDF with cut-off 6. 


[(Weber et al., 2019)]: https://pubs.rsc.org/en/content/articlehtml/2018/re/c7re00129k
[experiment1_strategic_molecules_random.py]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/3_screening/Screening_experiments/experiment1_strategic_molecules_random.py
[experiment1_random.py]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/3_screening/Screening_experiments/experiment1_random.py
[Index number (numbering within graph-tools)]: https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.Graph.vertex_index
[experiment2_strategic_molecuels_ranking.py]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/3_screening/Screening_experiments/experiment2_strategic_moelcules_ranking.py
[experiment3_DF_RPDF.py]: https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/3_screening/Screening_experiments/experiment3_DF_RPDF.py
[experiment3_RPDF.py]:https://github.com/Jana-Marie-Weber/strategic_molecules/blob/master/3_screening/Screening_experiments/experiment3_RPDF.py



