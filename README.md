Design of OpenMP-based Parallel Dynamic [PageRank algorithm] for measuring importance.

PageRank serves as an algorithm assessing the significance of nodes within a network through the assignment of numerical scores based on link structure. Its applications span web page ranking, identification of misinformation, traffic flow prediction, and protein target identification. The growing availability of extensive graph-based data has spurred interest in parallel algorithms for PageRank computation.

In real-world scenarios, graphs often evolve over time, with frequent edge insertions and deletions rendering recomputation of PageRank from scratch impractical, especially for small, rapid changes. Existing strategies optimize by iterating from the previous snapshot's ranks, reducing the required iterations for convergence. To further enhance efficiency, it becomes crucial to recalibrate the ranks of only those vertices likely to undergo changes. A common approach entails identifying reachable vertices from the updated graph regions and limiting processing to such vertices. However, if updates are randomly distributed, they may frequently fall within dense graph regions, necessitating the processing of a substantial portion of the graph.

To mitigate computational effort, one can incrementally expand the set of affected vertices from the updated graph region, rather than processing all reachable vertices from the initial iteration. Moreover, it is possible to skip processing a vertex's neighbors if the change in its rank is small and expected to have minimal impact on the ranks of neighboring vertices. Here, we introduce the Dynamic Frontier (DF) anf Dynamic Frontier with Pruning (DF-P) approaches for Updating PageRank, which addresses these considerations.

<br>


Below we plot the average time taken by Static, Naive-dynamic (ND), Dynamic Traversal (DT), our improved Dynamic Frontier (DF), and Dynamic Frontier with Pruning (DF-P) PageRank on 5 real-world dynamic graphs, with batch updates of size `10^-5|Eᴛ|` to `10^-3|Eᴛ|`. The labels indicate the speedup of each approach with respect to Static PageRank. DF PageRank is, on average, `8.0×`, `4.5×`, and `3.2×` faster than Static PageRank, and is, on average, `1.3×`, `1.1×`, and `1.5×` faster than DT PageRank, a widely used approach for updating PageRank on dynamic graphs. In contrast, DF-P PageRank is, on average, `26.2×`, `11.9×`, and `7.5×` faster than Static PageRank and, on average, `4.2×`, `2.8×`, and `3.6×` faster than DT PageRank on identical batch updates.

[![](https://i.imgur.com/YdjQWfH.png)][sheets-o1]

Next, we plot the Error comparison of Static, ND, DT, DF, and DF-P PageRank with respect to a Reference Static PageRank (with a tolerance `𝜏` of `10^−100` and limited to `500` iterations), using L1-norm.

[![](https://i.imgur.com/h2ZErIn.png)][sheets-o1]

Finally, we plot the strong scaling behaviour of DF and DF-P PageRank. With doubling of threads, DF PageRank exhibits an average performance scaling of `1.8×`, while DF-P PageRank exhibits an average performance scaling of `1.7×`.

[![](https://i.imgur.com/uahK7bg.png)][sheets-o2]


<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.


<br>
<br>


### Code structure

The code structure of Dynamic Frontier (DF), and Dynamic Frontier with Pruning (DF-P) PageRank is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/pagerank.hxx: PageRank algorithms
- inc/pagerankPrune.hxx: Dynamic Frontier with Pruning PageRank algorithms
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetricize.hxx: Graph Symmetricization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


