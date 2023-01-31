"""
"""

from utils import REVCOM

from typing import Tuple, Union

import sys
import os


class Sequence():
    """
    """

    def __init__(self, seqname, sequence) -> None:
        assert isinstance(seqname, str)
        assert isinstance(sequence, str)
        self._seqname = seqname
        self._sequence = sequence

    def __str__(self) -> str:
        return f">{self._seqname}\n{self._sequence}\n"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: sequance of {len(self)} characters>"

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index: Union[int, slice]) -> str:
        if isinstance(index, slice):
            return self.__class__(self._seqname, self._sequence[index])
        elif isinstance(index, int): 
            return self._sequence[index]
        else:
            raise TypeError(
                (
                    f"{self.__class__.__name__} subscripting requires {int.__name__} "
                    f"or {slice.__name__}, got {type(index).__name__}"
                )
            )
        
    def __setitem__(self, index: Union[int, slice], item: str) -> None:
        if not isinstance(item, str):
            raise TypeError(f"Expected {str.__name__}, got {type(item).__name__}")
        if isinstance(index, int):
            sequence = self._sequence[:index] + item + self._sequence[(index + 1):]
            self._sequence = sequence  
        elif isinstance(index, slice):
            start = 0 if index.start is None else index.start
            stop = len(self) if index.stop is None else index.stop
            sequence = self._sequence[:start] + item + self._sequence[stop:]

    def __iter__(self):
        return SequenceIterator(self)

    def _get_seqname(self) -> str:
        return self._seqname

    @property
    def seqname(self) -> str:
        return self._get_seqname()

    def _get_sequence(self) -> str:
        return self._sequence

    @property
    def sequence(self) -> str:
        return self._get_sequence()

    def _reverse_complement(self) -> str:
        return "".join([REVCOM[nt] for nt in self._sequence[::-1]])
    
    @property
    def reverse_complement(self) -> str:
        return self._reverse_complement()


class SequenceIterator():
    def __init__(self, sequence: Sequence) -> None:
        self._sequence = sequence
        self._i = 0

    def __next__(self):
        if self._i < len(self._sequence):
            e = self._sequence[self._i]
            self._i += 1
            return e
        raise StopIteration



