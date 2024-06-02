1.1 Boundary Initialization (David Baier)

With starting with the first boundary vertex, I iterate over all its neighbors and follow the ones, which are also a boundary edge. Then the edge lenght is computed by the distance of those two and stored separatly in an distintc vector. There is also an implementation with unordered_map, but the resulting circle is somehow messy. I guess this has to do with how the values are stored in the data structure. 
This is then repeated for all following boundary vertices. As soon as we hit an already available vertice, we stop the iteration and assume that all the vertices on the boundary are fetched. By computing the total length and then compute the ratio of the specific parts the circle is consturcted.

There where some problems with the unordered_map datastructure. But due to some time issues they were not resolved.

1.2 Iterative Solver (David Baier)
For the iterative solver, I run the number of times provided, each time computing the edge_weight. Iterating over each vertex and its neighbors the weights are computed and applied. After each vertex the summed weight and the summed weighted neigbor vertices are updated. 
Though the updates are made and the new positions are updated, all the vertices are removed. I assume this is to a updated position which results in a NaN. Due to time restrictions, i was not able to solve this problem.
