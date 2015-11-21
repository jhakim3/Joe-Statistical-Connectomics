Project Proposal - Joe and the Joes
===================
Question
-------------
We will consider, using cluster silhouetting, whether the distribution of the amount of clusters is more inherent to the species or the individuals. We will use the rats given, because there are three graphs available created using similar techniques.

>H0: The mean of the silhouette distributions between the individuals is not significantly different
Ha: The mean of the silhouette distributions between the individuals is significantly different

We will evaluate the difference in the silhouette profiles by calculating the RMSE between each curve, as well as comparing the means and standard deviations. We can acquire the latter by performing the clustering over several trials, giving us a distribution of optimal clusters as well as values for each cluster.

By establishing the similarity of the clustering, we should be able to apply this to the similarity of the overall architecture of the brain. 
Significance
--
Determining whether connectivity clustering is similar across species will provide a basis for further research into the implications of connectivity on organisms function and even possibly allow for the quantification of the effect of environmental differences on cognitive processes.
Innovation
-------------
In the world of graph theory and statistical connectomics, the quantification of the similarity between graphs has been an ongoing issue. Our experiment presents a method that evaluates similarity between graphs specific to rat brain networks using concepts in statistical analysis and clustering methods, namely the use of cluster silhouetting as a model.
Model
-------------
We will be modeling the different graphs as Stochastic Block Model random graphs. Doing this will not indicate the similarity of the clustering between the different individuals by itself, but it gives us an ability to interpret the similarity and difference. Namely, by indicating that the subgraphs behave as Erdős–Rényi model graphs, and that the off-diagonal blocks are more weakly associated than the diagonal ones, we can associate the similarity results with a similarity in connectivity of different brain regions.

Data
-------------
We will use the three rat graphs given in the Brain Architecture Management System Rat Connectome Project. This data set is the best available for comparing individuals from a population of individuals and compare properties between them.
Algorithm
-------------
We will model the brain connectivity using a stochastic block model where groups are determined by clustering (k means) to group our n nodes into k voronoi cells. While this is computationally intensive, it works well with large datasets because it is less expensive than many other algorithms. This would be even more of a determining factor if further investigation were to be done with the human brain, which is often represented with hundreds to thousands more nodes.

After grouping our nodes with the k-means method, we will use silhouetting to produce a graphical representation of how well the graph can be represented with some number of clusters. We will then compare the silhouettes using root mean squared error to establish the similarity or dissimilarity of the individual’s clustering tendencies. 
Expected Outcome
-------------
We expect to determine the level of similarity of the segmentation of the different regions of the brain, between individuals. This way, by establishing different levels of separations between different brain regions, the information about the distribution of clustering silhouette coefficients can be used as a sort of metric to demonstrate the levels of connectivity unique to an indivdiaul, or how an individual differs from the population.
Potential Pitfalls
-------------
We realize that there are certain problems in our model and algorithm of choice. The k means clustering algorithm exhibits intrinsic variability in its termination, and so given a particular number of clusters, there may be more than one classification of where the clusters are. Additionally, since we cluster based on number of edges, different groups could be incorrectly clustered together purely due to a similar number of edges. Additionally, we don’t use a stochastic block model for our data, as this proposal strays slightly from the pure topic of statistical connectomics.

In our model, the similarity comparison only uses one parameter - number of clusters, and so adding more parameters to the model would increase specificity. Additionally, the graphs may have an intrinsic similarity because they all map brains, and our model may not be robust to these similarities. Finally, the small sample size limits the statistical power of our results.