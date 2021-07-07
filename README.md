# Irregularity Index for Vector-valued Morphological Operators

Vector-valued morphological operators defined on total orders usually introduced irregularities because the topology of a total order is to weak to reproduce the topology induced by a metric. This repository contains the Julia's souce-code for measuring the irregularity issue using the relative gap between the generalized sum of pixel-wise distances and the Wasserstein metric. Currently, the source-code allows only computing the irregularity of color images using the Euclidean distance and with the parameters p=1 or p=2.

## Getting Started

This repository contain the Julia's source-codes for quantifying the irregularity issue of vector-valued morphological operators. The Jupyter-notebook of the computational experimens are also available in this repository.

## Usage

First of all, call the irregularity tools module using:

```
include("IrregularityTools.jl")
```

The modole uses the following libraries: Images, ImageMorphology, Statistics, StatsBase, Distances, LinearAlgebra, ProgressMeter, OptimalTransport, Tulip, JuMP, and Clp. Add these libraries before calling the irregularity tools module.

### Global Irregularity Index:

The global irregularity index evaluate the gap between the generalized sum of pixel-wise distances and the Wasserstein metric computed using all the pixels of an image. Due to the computational cost, the global irregularity index can be computed only for tiny images. The global irregularity index can be computed as follows:
```
global_irregularity(I,J)
```
where I is the input image and J is the output image of a morphological operator. The previous command computes the irregularity index by solving the optimal transport problem using a linear programming solver (JuMP + Clp). Alternatively, the global irregularity index can be computed using the Sinkhorn method or its stabilized version, both available at OptimalTransport library.
```
global_irregularity(I,J,"sinkhorn",1.e-2)
```
or
```
global_irregularity(I,J,"stabilized_sinkhorn",1.e-3)
```

### Local Irregularity Index:

The local irregularity index is obtained by computing the Wasserstein metric in local windows of size $W \times W$ and aggregating the results in a single value. The local irregularity index provides a lower bound for the global irregularity index. The local irregularity index with local windows of size $16 \times 16$ can be computed as follows:
```
local_irregularity(I,J)
```
where I is the input image and J is the output image of a morphological operator. The size of the local windows can be given as an additional parameter:
```
local_irregularity(I,J,W)
```
The previous commands compute the irregularity index by solving the optimal transport problem using the stabilized Sinkhorn method from available at OptimalTransport library. The Sinkhorn method can also be used as follows:
```
local_irregularity(I,J,W,"sinkhorn",1.e-2)
```
Alternatively, the optimal transport can be solved analytically using JuMP + Clp as follows:
```
local_irregularity(I,J,W,"JuMP")
```

### See examples in the Jupyter notebook files!
