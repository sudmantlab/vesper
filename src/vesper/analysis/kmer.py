from __future__ import annotations

from typing import Set, Dict
# from ..models.variants import Variant, VariantAnalysis, SVType
# from ..models.reads import AlignedRead


# this is better done with jellyfish/khmer, was only for proof of concept

class KmerDecomposition:
    """Represents a sequence decomposed into kmers."""
    k: int
    kmers: Set[str]
    kmer_counts: Dict[str, int]
    kmer_frequencies: Dict[str, float]

    def __init__(self, sequence: str, k: int):
        self.k = k
        self._decompose(sequence)

    def __repr__(self):
        """Return a string representation of the kmer decomposition."""
        return f"KmerDecomposition(k={self.k}, n={len(self.kmers)}), kmers = {self.kmer_frequencies}"
    
    def __len__(self):
        """Return the number of kmers in the decomposition."""
        return len(self.kmers)

    def _decompose(self, sequence):
        """Decompose the sequence into kmers."""
        kmers = []
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            kmers.append(kmer)
        self.kmers = set(kmers)
        self.kmer_counts = {kmer: kmers.count(kmer) for kmer in self.kmers}
        self.kmer_frequencies = dict(sorted(
            {kmer: round(self.kmer_counts[kmer] / len(kmers), 5) for kmer in self.kmers}.items(),
            key=lambda x: x[1],
            reverse=True
        ))

    @classmethod
    def difference(cls, k1: 'KmerDecomposition | set[str]', k2: 'KmerDecomposition | set[str]'):
        """Returns the set of kmers not shared between two sets."""
        set1 = k1.kmers if isinstance(k1, KmerDecomposition) else k1
        set2 = k2.kmers if isinstance(k2, KmerDecomposition) else k2
        return set1 - set2
    
    @classmethod
    def intersection(cls, k1: 'KmerDecomposition | set[str]', k2: 'KmerDecomposition | set[str]'):
        """Returns the set of kmers shared between two sets."""
        set1 = k1.kmers if isinstance(k1, KmerDecomposition) else k1
        set2 = k2.kmers if isinstance(k2, KmerDecomposition) else k2
        return set1 & set2
    
    @classmethod
    def union(cls, k1: 'KmerDecomposition | set[str]', k2: 'KmerDecomposition | set[str]'):
        """Returns the set of kmers in either set."""
        set1 = k1.kmers if isinstance(k1, KmerDecomposition) else k1
        set2 = k2.kmers if isinstance(k2, KmerDecomposition) else k2
        return set1 | set2
    
    @classmethod
    def jaccard_index(cls, k1: 'KmerDecomposition | set[str]', k2: 'KmerDecomposition | set[str]'):
        """Returns the Jaccard index of two sets."""
        set1 = k1.kmers if isinstance(k1, KmerDecomposition) else k1
        set2 = k2.kmers if isinstance(k2, KmerDecomposition) else k2
        intersection = set1 & set2
        union = set1 | set2
        return len(intersection) / len(union)
    
    @classmethod
    def cosine_similarity(cls, k1: 'KmerDecomposition | set[str]', k2: 'KmerDecomposition | set[str]'):
        """Returns the cosine similarity of two sets."""
        set1 = k1.kmers if isinstance(k1, KmerDecomposition) else k1
        set2 = k2.kmers if isinstance(k2, KmerDecomposition) else k2
        intersection = set1 & set2
        return len(intersection) / (len(set1) * len(set2))**0.5
