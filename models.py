from fm_index_query import FMIndexQuery     
from global_alignment_algo import globalAlignment
from local_alignment_algo import localAlignment


class DNASequence:
    def __init__(self, seq: str, ID: str, description: str = ""):
        self._seq = seq
        self._id = ID
        self._description = description

    def get_subsequence(self, start: int, end: int):
        if start < 0 or end > len(self._seq):
            raise ValueError("Subsequence indices out of range")
        return self._seq[start:end]

    def get_GC_content(self):
        gc_count = self._seq.count("G") + self._seq.count("C")
        return (gc_count / len(self._seq)) * 100

    def get_length(self):
        return len(self._seq)

    @property
    def seq(self):
        return self._seq

    @property
    def id(self):
        return self._id

    @property
    def description(self):
        return self._description


class MitochondrialDNA(DNASequence):
    def __init__(self, seq: str, ID: str, description: str = ""):
        super().__init__(seq, ID, description)


class MotifFinder:
    def __init__(self, motif_seq: str):
        self._motif_seq = motif_seq

    def count_occurrences(self, target_seq: str):
        return FMIndexQuery(target_seq, self._motif_seq)[0]

    def search_motif(self, target_seq: str):
        return FMIndexQuery(target_seq, self._motif_seq)[1]

class SequenceAlignment:
    def __init__(self, seq1: str, seq2: str):
        self._seq1 = seq1
        self._seq2 = seq2

    def align_sequences(self, gap_pen=-2, match=1, mismatch=-1, algo: str = "global"):
        if algo == "global":
            seq1_gapped, comparison, seq2_gapped = globalAlignment(
                self._seq1, self._seq2, gap_pen, match, mismatch
            )[0]
        elif algo == "local":
            seq1_gapped, comparison, seq2_gapped = localAlignment(
                self._seq1, self._seq2, gap_pen, match, mismatch
            )[0]
        else:
            raise ValueError("Unknown alignment algorithm.")
        return seq1_gapped, comparison, seq2_gapped

    def get_alignment_scores(self, gap_pen=-2, match=1, mismatch=-1, algo: str = "global"):
        if algo == "global":
            score = globalAlignment(self._seq1, self._seq2, gap_pen, match, mismatch)[1]
        elif algo == "local":
            score = localAlignment(self._seq1, self._seq2, gap_pen, match, mismatch)[1]
        else:
            raise ValueError("Unknown alignment algorithm.")
        return score
