# Dictionary Learning with Rapid Orthogonal Matching Pursuit

## Short Description
Fast OMP sparse representation processing by utilizing matrix operations, which are quicker than loops.

## Publication
DOI: [10.1109/igarss39084.2020.9323532](https://doi.org/10.1109/igarss39084.2020.9323532)

## Paper Abstract
Orthogonal Matching Pursuit (OMP) has proven itself to be a significant algorithm in image and signal processing domain in the last decade to estimate sparse representations in dictionary learning. Over the years, efforts to speed up the OMP algorithm for the same accuracy has been through variants like generalized OMP (g-OMP) and fast OMP (f-OMP). All of these algorithms solve OMP recursively for each signal sample among S number of samples. The proposed rapid OMP (r-OMP) runs the loop for N atoms, simultaneously estimating for all samples, and, in a real scene since N ≪ S, the proposed approach speeds up OMP by several orders of magnitude. Experiment on a real scene with a popular dictionary learning algorithm, K-SVD, show that the proposed r-OMP completes K-SVD in ≈4% of the computational time compared to using OMP.
