# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""A weighted event index with logarithmic updates and selection (Fenwick tree)."""

import math
from functools import lru_cache

import numpy as np
from numpy.typing import NDArray


@lru_cache(maxsize=8)
def _update_paths(size: int):
    """Cache Fenwick ancestors in the same leaf-by-leaf order as update().

    Bounded caching shares integer paths between equal-sized trees, without
    retaining per-run rates or random state. Arrays are private and read-only.
    """
    paths = np.zeros((size, max(1, size.bit_length())), dtype=np.intp)
    indices = np.arange(1, size + 1, dtype=np.intp)
    for column in range(paths.shape[1]):
        valid = indices <= size
        paths[valid, column] = indices[valid]
        indices = indices + (indices & -indices)
    counts = np.count_nonzero(paths, axis=1)
    flattened = paths[paths > 0]
    for array in (paths, counts, flattened):
        array.flags.writeable = False
    return paths, counts, flattened


class RateTree:
    """Store non-negative event rates without scanning every channel per event.

    Each tree entry stores a partial sum. Updating or selecting a channel visits
    O(log n) entries. This changes bookkeeping cost, not Gillespie probabilities.
    """

    def __init__(self, values: NDArray[np.float64]) -> None:
        self.values = np.asarray(values, dtype=float).copy()
        if np.any(~np.isfinite(self.values)) or np.any(self.values < 0):
            raise ValueError("Event rates must be finite and non-negative.")
        self.tree = np.zeros(len(self.values) + 1)
        self.tree[1:] = self.values
        self.positive = int(np.count_nonzero(self.values))
        for index in range(1, len(self.tree)):
            parent = index + (index & -index)
            if parent < len(self.tree):
                self.tree[parent] += self.tree[index]

    @classmethod
    def from_sequential_updates(cls, values: NDArray[np.float64], order=None):
        """Reproduce zero-tree sequential insertion without Python update loops.

        A bulk constructor can sum floating-point values in a different order.
        Here unbuffered add.at retains the old per-leaf ancestor-addition order,
        including the site-then-polymer insertion order of the chemical tree.
        This preserves the previous seeded selection/time arithmetic exactly.
        """
        result = cls.__new__(cls)
        result.values = np.asarray(values, dtype=float).copy()
        if (result.values.ndim != 1 or np.any(~np.isfinite(result.values))
                or np.any(result.values < 0)):
            raise ValueError("Event rates must be a finite non-negative vector.")
        size = len(result.values)
        paths, counts, indices = _update_paths(size)
        rates = result.values
        if order is not None:
            order = np.asarray(order)
            if (order.shape != (size,) or order.dtype.kind not in "iu"
                    or not np.array_equal(np.sort(order), np.arange(size))):
                raise ValueError("Insertion order must be a permutation of leaf indices.")
            ordered_paths = paths[order]
            indices = ordered_paths[ordered_paths > 0]
            rates, counts = rates[order], counts[order]
        result.tree = np.zeros(size + 1)
        np.add.at(result.tree, indices, np.repeat(rates, counts))
        result.positive = int(np.count_nonzero(result.values))
        return result

    @property
    def total(self) -> float:
        """Return the root sum, treating an entirely disabled family as zero."""
        if self.positive == 0:
            return 0.0
        total, index = 0.0, len(self.values)
        while index:
            total += float(self.tree[index])
            index -= index & -index
        return total

    def update(self, index: int, value: float) -> None:
        """Replace one rate and update only the partial sums containing it."""
        if not math.isfinite(value) or value < 0:
            raise ValueError("Event rates must be finite and non-negative.")
        previous = float(self.values[index])
        if value == previous:
            # Preserve the supplied leaf value (including signed zero), while
            # avoiding ancestor walks for unchanged local rates.
            self.values[index] = value
            return
        self.positive += int(value > 0) - int(previous > 0)
        self.values[index] = value
        difference = value - previous
        index += 1
        while index < len(self.tree):
            self.tree[index] += difference
            index += index & -index

    def select(self, weight: float) -> int:
        """Locate the channel containing a draw in [0, total)."""
        if not 0 <= weight < self.total:
            raise ValueError("Selection weight must lie below a positive total rate.")
        index = 0
        bit = 1 << (len(self.values).bit_length() - 1)
        while bit:
            candidate = index + bit
            if candidate < len(self.tree) and self.tree[candidate] <= weight:
                weight -= float(self.tree[candidate])
                index = candidate
            bit >>= 1
        if index >= len(self.values) or self.values[index] <= 0:
            raise ArithmeticError("Weighted index drifted; a cache rebuild is required.")
        return index
