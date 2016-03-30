# CTSTA
Characterisation of microstructures


1. The source code is released accompanying the publication of the paper “Jie Liu, Gerald G. Pereira, Qingbin Liu, Klaus Regenauer-Lieb, 2016, Computational challenges in the analyses of petrophysics using microtomography and upscaling: A review, Computers & Geosciences, 89, 107–117, doi:/10.1016/j.cageo.2016.01.014.”
Please cite this paper if you use the tool to do characterisation in your research.


2. The code calls a function of mathematical library of Intel, thus it is preferred to compile it using Intel-compiler and use “-mkl” to call the math-lib.


3. Some basic concepts

3.1 Connection: in a binary dataset after segmentation, two voxels with a common plane and of the same material is connected; two voxels with a common line and of the same material is not connected.

3.2 Cluster: a cluster is a group of voxels connected one by one and of the same material. A cluster can be very large and complex; also can be very small and simple – an isolated voxel with different material is defined as a cluster.

3.3 Labelling clusters: is a process of giving all voxels within the same cluster the same label. The largest cluster (with maximum voxel number) is cluster 1, the second largest one is cluster 2, and so on.

3.4 Percolating: We define a direction is percolating when one boundary of the parallelepiped in the direction has the same cluster label with its opposite boundary. This means the two end boundaries of the direction are connected by the same cluster and the model is permeable.

