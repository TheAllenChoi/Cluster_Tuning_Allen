## Overview

This repo is for collecting experiments in how to use resampling to tune hyperparameters in **unsupervised** learning.

Unlike in *supervised* learning, in unsupervised learning we cannot identify a singular metric that measures the fit of the model. Instead, we seek to find principles of "good fit" that can be compared across hyperparameter values.

Possible resampling structures for clusting include:

### Cluster consistency across reruns with random conditions

Some clustering algorithms, such as $k$-means, include randomized starting conditions. Can we exploit this as a form of "resampling", and look for consistency of certain cluster properties across runs?

(`initial_conditions` folder)

### Cross-validation and classic metrics

Many common clustering metrics are guaranteed to shrink as $k$ increases; for example, the within-sum-of-squares metric is smaller when $k$ is larger. For this reason, a classic approach is to use the "elbow method" to identify when the "gains" start to diminish. We believe this approach is ambiguous and unsatisfying.

However, it may be interesting to see how these metrics behave over cross-validation folds. Perhaps they are more variable, or more stable, across folds at the correct $k$? Perhaps a metric can be created to compare each fold to a "correct" answer?

(`cross_validation` folder)

### Cluster consistency across subsamples

I believe this is the most promising avenue.

In principle, if a clustering structure is "real", it should emerge in any random sample of data. Can we use some kind of bootstrapping or subsampling approach on our data to see if the same structure emerges in different samples?

The challenge here is measuring "same structure". Do we use cluster centers? Cluster membership (but then how to compare across subsamples with different items in them)? Other cluster properties?

(`subsampling` folder)

## Phase 1

In Phase 1, we focus only on **number of clusters** as the hyperparameter of interest, and we restrict our experiments to **k-means clustering** and **hierarchical clustering**.
