# Weisfeiler-Lehman Subtree Kernel

This programme computes the [Weisfeler-Lehman Subtree Kernel](https://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf) on two undirected, weighted graphs. It can also be used to compute embeddings of two graphs produced by the Weisfeiler-Lehman test of graph isomorphism.

## Third-Party

The label compression is based on the MD5 hashing from [Brumme's Hashing Library](https://create.stephan-brumme.com/hash-library/).

## Installation

On a bash command, at the root of the folder:

```bash
cd WLKernel
make
```

## Usage

In folder WLKernel:
```
./wlkernel <File_1> <File_2> [options]
```
Computes the similarity between the graphs in `<File_1>` and `<File_2>`, based on the normalised dot product of their embeddings produced by a Weisfeiler-Lehman test.
* If extension of `<File_i>` is `emb`, the programme assumes that the file contains an already computed embedding.

<img align="right" src="https://github.com/acaen/MARGOT/blob/master/Pics/toyGraph.png" width="15%" height="15%">

* Otherwise, it computes the graph embedding using a Weisfeiler-Lehman test of depth 2. Initial node features are their weighted degrees. Graphs must be in a weighted edgelist format, with integer nodes (see Figure).

Options:
* `-d <k>` : To change the depth of the Weisfeiler-Lehman test to `<k>`.

:warning: This has no impact if both graph embeddings are read directly.

* `-save` : If the embedding of one or both graph(s) are built, the programme saves it in a file `<File_i>.emb`.
* `-nonorm`: To not normalise the dot product of the two graph embeddings (otherwise it is normalised by the max of both norms).
* `-feat_k <FileToFeat>`: Path to a file containing node features to use instead of node degrees. Must contain as many rows as there are nodes in the format `node feat`. Features must be integers.

:warning: This has no impact if the corresponding graph embedding is read directly.

## Example

```
./wlkernel ../Data/g1.txt ../Data/g2.txt -d 1 -feat1 ../Data/feat1.txt -feat2 ../Data/feat2.txt -nonorm
```

Graphs `g1` and `g2` with respective node features `feat1` and `feat2` are graphs from Fig. 2 of [Shervashidze's paper](https://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf). With a Wesfeiler-Lehman test of depth 1 and no normalisation, their kernel is found to be 11.
