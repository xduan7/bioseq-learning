"""
File Name:          genome_sequence.py
Project:            bioseq-learning

File Description:

"""
import logging
from enum import Enum

from typing import Any, Dict,  Optional, Sequence, Tuple

from Bio import Seq


_LOGGER = logging.getLogger(__name__)


class AnnotationType(Enum):
    CodingRegion = 'CodingRegion'
    ConservedDomain = 'ConservedDomain'


class _Annotation:
    def __init__(
            self,
            range: Tuple[int, int],
            attributes: Dict[str, Any],
            annotation_type: AnnotationType,
    ):
        self._start: int = range[0]
        self._end: int = range[1]
        self._attributes: Dict[str, Any] = attributes
        self._annotation_type: AnnotationType = annotation_type
        # other common attributes for genome sequence annotations


class CodingRegion(_Annotation):
    def __init__(
            self,
            range: Tuple[int, int],
            attributes: Dict[str, Any],
    ):
        super().__init__(range, attributes, AnnotationType.CodingRegion)
        # other attributes specific to coding region


class ConservedDomain(_Annotation):
    def __init__(
            self,
            range: Tuple[int, int],
            attributes: Dict[str, Any],
    ):
        super().__init__(range, attributes, AnnotationType.ConservedDomain)
        # other attributes specific to conserved domain (e.g. hit types)


class GenomeContigSeq(Seq):

    def __init__(
            self,
            genome_contig_seq_path: str,
            genome_contig_cr_tsv_path: Optional[str],
            genome_contig_cd_csv_path: Optional[str],

    ):
        # TODO: assert that there is only one sequence for the file?
        self.coding_regions: Sequence[_Annotation]
        self.conserved_domains: Sequence[_Annotation]
