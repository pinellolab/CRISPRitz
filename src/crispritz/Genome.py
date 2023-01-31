"""
"""

from Sequence import Sequence

from typing import List, Union

import os


class FastaReader():
    """
    """

    def __init__(self, fname) -> None:
        assert os.path.isfile(fname)
        self._fname = fname

    def read(self) -> Union[List[Sequence], Sequence]:
        # skip comment lines
        handle = open(self._fname, mode="r")
        while True:
            line = handle.readline().strip()
            if not line:
                return  # empty file?
            if line.startswith(">"):
                break  # data start here
        sequences = []
        while True:
            if not line.startswith(">"):
                raise ValueError(f"FASTA sequence name should be preceded by '>'")
            seqname = line[1:].strip().split()[0]
            lines = []
            line = handle.readline().strip()
            while True:
                if not line:
                    break  # empty sequence ?
                if line.startswith(">"):
                    break  # sequence end
                lines.append(
                    line.strip().replace(" ", "").replace("\r", "").replace("\n", "")
                )
                line = handle.readline().strip()
            sequences.append(Sequence(seqname, "".join(lines)))
            if not line:
                break  # EOF
        if len(sequences) == 1:
            return sequences[0]
        return sequences


class Genome():
    """
    """

    def __init__(self, sequences: List[Sequence]) -> None:
        self._sequences = {s.seqname: s for s in sequences}
        self._sequences_list = list(self._sequences.values())

    def __getitem__(self, index) -> Sequence:
        if index not in self._sequences.keys():
            raise KeyError(f"Forbidden sequence name ({index})")
        return self._sequences[index]

    def __len__(self) -> int:
        return len(self._sequences_list)

    def __iter__(self):
        return GenomeIterator(self)


class GenomeIterator():
    def __init__(self, genome) -> None:
        self._genome = genome
        self._i = 0

    def __next__(self) -> Sequence:
        if self._i < len(self._genome):
            e = self._genome._sequences_list[self._i]
            self._i += 1
            return e
        raise StopIteration
