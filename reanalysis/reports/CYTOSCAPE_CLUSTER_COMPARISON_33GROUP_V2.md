# Cytoscape clustering comparison for the fixed 33-group Kla–DDR analysis

## Scope

- The analysis is restricted to the fixed 33 paired-group project scope and the
  507 proteins in Kla ∩ DDR.
- Protein identity is keyed by isoform-stripped UniProt `BaseAccession`.
  `GeneSymbol` is retained only as a display label.
- The STRING v12 network was queried at a confidence cutoff of 0.70 with zero
  additional interactors. STRING mapped 506 proteins; `A8MVJ9` was retained as
  an isolated node because it has no STRING mapping.
- The final network contains 507 nodes and 4,458 edges.

## Three versions retained

1. `Kla_DDR_507_STRING_MCL_i3_shared_STRING_topology`
   - Network-driven reference version.
   - clusterMaker2 MCL, inflation = 3.0, STRING combined score as the edge
     weight.
   - 65 modules, 446 clustered proteins and 61 unclustered proteins.
   - This is the version that does **not** use the later manually defined
     pathway-clustering request.
2. `Kla_DDR_507_pathway_score_k7_seed25_shared_STRING_topology`
   - Hypothesis-driven seven-dimensional version.
   - Features: HR, NHEJ, AEJ, BER, NER, MMR and FA.
   - 189 proteins have an all-zero seven-pathway score.
3. `Kla_DDR_507_pathway_score_k9_seed25_shared_STRING_topology`
   - Hypothesis-driven nine-dimensional version.
   - Adds `Chromatin interaction` and `Other support`; these are auxiliary
     categories rather than canonical DNA-repair pathways.
   - 22 proteins have an all-zero nine-pathway score.

For the two hypothesis-driven versions, each protein was encoded as −1
(suppressing), 0 (unassigned) or +1 (promoting). K-means was run in R with
`set.seed(25)`, `nstart = 250`, `iter.max = 1000` and requested `k = 7` or
`k = 9`. The resulting clusters are pathway-score **patterns**. They must not be
described as seven or nine mutually exclusive biological pathways.

The seven-dimensional solution has a mean silhouette width of 0.605, compared
with 0.400 for the nine-dimensional solution. Thus the seven-dimensional
solution is more clearly separated statistically, whereas the nine-dimensional
solution assigns substantially more proteins to at least one category.

## Shared topology layout

The three primary views use exactly the same per-protein X/Y coordinates. A
force-directed layout is calculated once from the 4,458 STRING edges and then
copied by `BaseAccession` to the MCL, seven-pattern and nine-pattern views.
Therefore, switching between the three networks changes node color and the
all-zero marker shape, but not protein position.

The graph has one 446-node connected component, four two-node components and 53
singletons. The 446-node component retains its force-directed geometry. The
smaller disconnected components are packed into a compact shelf below it so
that isolated nodes do not create a mostly empty canvas. This only changes the
placement of components that have no edge path between them; it does not invent
edges or alter topology within the main connected component.

The earlier `Group Attributes Layout` views are retained only as
`*_grouped_circle_layout` comparisons. They intentionally arrange clusters as
circles and must not be used as the topology-preserving primary result.

## Files

- Cytoscape session:
  `reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1/kla_ddr_507_cluster_comparison_shared_string_topology_33groups_v2.cys`
- Node import table:
  `reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1/cytoscape_node_import_table.csv`
- Normalized STRING nodes:
  `reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1/cytoscape_string_nodes_507.csv`
- Normalized STRING edges:
  `reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1/cytoscape_string_edges_4458.csv`
- Cluster assignments:
  `reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1/pathway_cluster_assignments_seed25.csv`
- Cluster centers and diagnostics are stored in the same table directory.
- PNG, SVG and PDF exports are stored in:
  `reanalysis/results/figures/kla_ddr_cytoscape_pathway_clusters_33groups_v1/`
- Primary PNG files are exported at 8,000 px width. SVG and PDF are vector
  formats and should be used when unlimited zoom or manuscript editing is
  needed.

## Fastest way to inspect the result

1. Open Cytoscape 3.10.4.
2. Choose `File > Open Session` and open
   `kla_ddr_507_cluster_comparison_shared_string_topology_33groups_v2.cys`.
3. In the Network panel, select one of the three networks listed above.
4. In the Node Table, inspect:
   - `MCL_cluster_i3`
   - `Score_HR` through `Score_Other_support`
   - `K7_Cluster`, `K7_ClusterLabel`
   - `K9_Cluster`, `K9_ClusterLabel`
5. In Style, labels are available from `DisplayLabel` but are hidden in the
   overview to avoid showing 507 overlapping labels.

## Reproduce the seven- or nine-pattern version in the Cytoscape GUI

The exact published labels should be calculated in R and imported, because
clusterMaker2's K-means GUI does not expose a random-seed setting.

1. Run:

   ```bash
   Rscript reanalysis/scripts/cluster_kla_ddr_pathway_scores_33groups.R
   ```

2. In Cytoscape, select the 507-node network.
3. Choose `File > Import > Table from File`.
4. Select `cytoscape_node_import_table.csv`.
5. Import the table into the selected network and match:
   - file key: `BaseAccession`
   - network key: `BaseAccession`
6. For the seven-dimensional version, set Style mappings:
   - Fill Color: discrete mapping from `K7_cluster_label`
   - Shape: discrete mapping from `K7_assignment_status`
   - `Assigned` = ellipse; `All-zero` = diamond
7. For the nine-dimensional version, use the corresponding `K9_...` columns.
8. Do not recalculate a layout for the seven- or nine-pattern primary view.
   Copy the topology reference coordinates by `BaseAccession`, for example with
   `Layout > Copycat Layout`:
   - source and target key: `BaseAccession`
   - source network:
     `Kla_DDR_507_STRING_MCL_i3_shared_STRING_topology`
   - target network: the seven- or nine-pattern network
9. Use `cytoscape_cluster_color_key.csv` to reproduce the saved colors.

`Layout > Group Attributes Layout` may be used only to reproduce the optional
circle comparison. It removes the STRING topology coordinates and should not be
used for the primary network figure.

## Manual save and export in Cytoscape

1. After checking the three primary views, choose `File > Save As` and save the
   complete Cytoscape session as a `.cys` file. Use `Save As`, rather than
   `Save`, when a new versioned filename is required.
2. Before exporting each network, select it in the Network panel and choose
   `View > Fit Content`.
3. Choose `File > Export > Network View as Graphics`.
4. Use SVG or PDF for a manuscript or unlimited zoom. For a raster copy, choose
   PNG and set a large export width; the scripted primary PNG files use 8,000
   px.
5. Keep node labels hidden in the full 507-node overview. To inspect individual
   proteins interactively, zoom in and enable the `DisplayLabel` node-label
   mapping temporarily.
6. If a pathway-score view is accidentally rearranged, do not rerun a
   force-directed or grouped-circle layout on it. Reapply `Copycat Layout` from
   the shared STRING-topology network using `BaseAccession` as both keys.

If only an exploratory GUI result is needed, clusterMaker2 can run K-means
directly:

1. Open `Apps > clusterMaker Cluster Attributes > K-Means Cluster`.
2. For the seven-dimensional run, choose the seven `Score_...` node columns and
   set the number of clusters to 7.
3. For the nine-dimensional run, choose all nine `Score_...` columns and set the
   number of clusters to 9.
4. Do not use the GUI-only result as the exact reproducible result because its
   random initialization is not fixed to seed 25.

## Reproduce the topology-only MCL reference

1. Import the 507 BaseAccession values with stringApp using:
   - species: *Homo sapiens*
   - confidence cutoff: 0.70
   - additional interactors: 0
2. STRING should recognize 506 values. Retain `A8MVJ9` as a separate isolated
   node.
3. Open `Apps > clusterMaker Cluster Network > MCL Cluster`.
4. Select STRING combined score as the edge weight.
5. Set inflation to 3.0 and save the result to a node column such as
   `MCL_cluster_i3`.
6. Apply a topology-preserving force-directed layout and color nodes by
   `MCL_cluster_i3`; leave unclustered proteins gray.

Changing MCL inflation does not set an exact cluster count. In the tested
network, inflation 2.5, 3.0, 3.5 and 4.0 produced 46, 65, 77 and 78 modules,
respectively. Therefore, a request for exactly seven or nine clusters must use
the pathway-score matrix rather than changing the MCL inflation value.
