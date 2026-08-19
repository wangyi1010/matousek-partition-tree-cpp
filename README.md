# Matoušek Partition Tree

## 1. The problem

Given a fixed set of points in the plane, preprocess them once so that each new
halfplane query can quickly answer one question:

> How many input points lie inside this halfplane?

A brute-force query checks every point. A partition tree instead counts whole
groups when possible and examines only the groups crossed by the query
boundary.

## 2. The algorithm

1. Convert the input points into lines using point-line duality.
2. Build a verified weighted cutting and derive a finite test set.
3. Use exponential reweighting to extract small point groups contained in
   triangles.
4. Recurse on those groups to build the partition tree.
5. For each query, count groups fully inside, skip groups fully outside, and
   recurse only into crossed groups.

## 3. Complexity

For any fixed positive approximation parameter, the theorem gives

```math
T_{\mathrm{query}}(n)
=
O\!\left(n^{1/2+\varepsilon}\right)
```

and

```math
T_{\mathrm{build}}(n)
=
O(n\log n).
```

[Read the complete mathematical derivation →](docs/two_dimensional_partition_theorem_math.md)

## 4. Visual example

![Brute-force scanning compared with partition-tree querying](assets/problem_comparison.svg)
