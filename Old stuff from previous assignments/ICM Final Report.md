***Introduction to Computational Medicine - Final Report***

Joe and the Joes

**Abstract:**

In this report, we will pose the question of whether the three rat brains in the Brain Architecture Management System Rat Connectome Project database are similar enough in terms of associativity to be described as similar brains. To do this, we intend to iterate through different cluster-sized stochastic block models (SBM), and at each step, determine if any of the three graphs are sampled from that sized SBM. If the three graphs all conform to the same set of SBMs, we will conclude that they are of similar associativity.

**Question:**

Informal: 

The problem we're trying to solve is whether all three rat brains can be modeled by the same size (number of clusters) stochastic block model, or whether they cannot. This will indicate the overall similarity in the brain associativity.

Formal: 

Given Gi, i ∈ {1,2,3} where Gi ∈Gis the space of all graphs that are directed and loopy and have 503 vertices and unweighted edges, and Gi denotes the adjacency graphs for three rat brains,
H0: Gi ~ SBMk
Ha: Gi ~ SBMj j ≠ k

Where k,h ∈ {1…5} denote number of blocks in the stochastic block model (SBM).

Test Statistic:
 ![](http://i.imgur.com/JdKoTI4.png)

Where Λ is a maximum likelihood ratio between the null and alternative hypotheses.

We will generate the null and alternative distributions using a Monte Carlo simulation on the respective SBM block model sizes.

For each k, we will produce a result for a given p < 0.05 whether each graph is distributed as a size k SBM. Each graph will have a set of k sized SBMs where they match, and we will match those sets of k between graphs. If there are elements in common, we will conclude that the rat brains are of the same associativity.

**Pseudocode:**

Inputs: rat brain files. Outputs: singificance levels for each combination of cluster sizes
Extract rat matrices from database
For each rat, create an adjacency matrix for their graphs
For clusters of size n and m,
Use spectral clustering algorithm to create two vertex sets for each rat, sized n and m
Find maximum likelihoods for each edge in each cluster
Use to create a random matrix sampled from SBM of each size
Compute the likelihood ratio of these random graphs being sampled from SBMs of each size
	Apply heuristic penalty for larger block size: subtract ~ number of blocks^2
Evaluate when likelihood ratio of rat graphs is greater or less than this ratio
	If all rats are distributed from one sized SBMs statistically, conclude they are similarly associative.

Spectral Clustering Pseudocode
Input: the adjacency matrix (W), the number of clusters (k), Type of normalization (Type)
Output: Clustercell
First iterate through the adjacency matrix to convert any element greater than or equal to 1 to 1
Then calculate the degree matrix from the sum of the simplified matrix
From the degree matrix and the simplified matrix, compute the normalized Laplacian
If specified, normalization Laplacian using algorithms specified in Shi and Malik (2000) (Type = 2) or Jordan and Weiss (2002) (Type = 3)
Compute the eigenvectors based off of the k smallest eigenvalues
Use k-means to cluster the eigenvectors
Place each cluster in a separate array and place all arrays into Clustercell

**Data:**

We will use the three rat graphs given in the Brain Architecture Management System Rat Connectome Project. This data set is the best available for comparing individuals from a population of individuals and compare properties between them.

Model

We will be modeling the different graphs as Stochastic Block Model random graphs. In order to do this, given the fact that the rat brains are weighted random graphs, we would need to apply some thresholding to make each edge’s connection binary instead of an integer. Hypothetically, this is not optimal for biological data such as brain connectomes, as other non-bernoulli type connection schemes generally are better models (scale-free networks, etc.) However, we chose to use SBMs for simplicity.


**Results:**

![](http://i.imgur.com/FOEPDx6.jpg)
![](http://i.imgur.com/JFu1FIh.jpg)	
![](http://i.imgur.com/w9TrBvd.jpg)

Rats 1,2,3: Significance levels for difference in block sizes. Since no size (1-5) is consistently the best fit for each rat, there is no overlap and therefore we cannot conclude that the rats are of similar associativity. For instance, rat #1 had H0 cluster size of 2 and significantly not 3, 4, or 5 (using the figure, and disregarding cluster size 1), however this was not observed for rats #2 or #3, (with rats #2 and #3, Ha also worked) etc. 


**Discussion:**

Although our sample size is low, if we take our results to indicative of a larger phenomena, clustering is not similar across species. This fact provokes further study into the significance of these clustering differences and how they translate to individual rat development differences. However, the conclusion that the rat brains are not similar is not accessible, since we only compared them using a very specific metric. In addition, we made several assumptions about the distribtion of the graph itself, including the statement that all the edge probabilities were independent, and identically distributed within a cluster and between clusters. However, in biological settings, independence is rare and perhaps a more scale-free graph structure would be a more accurate description. 

**Significance:**

Determining whether connectivity clustering is similar across species will provide a basis for further research into the implications of connectivity on organisms function and even possibly allow for the quantification of the effect of environmental differences on cognitive processes.

Sources:
Shi and Malek (2000) http://www.cs.berkeley.edu/~malik/papers/SM-ncut.pdf
Jordan and Weiss (2002) http://ai.stanford.edu/~ang/papers/nips01-spectral.pdf
Spectral Clustering Code (Ingo Buerk 2011/2012)
 http://www.mathworks.com/matlabcentral/fileexchange/34412-fast-and-efficient-spectral-clustering

