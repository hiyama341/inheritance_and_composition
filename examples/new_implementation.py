from dataclasses import dataclass
from math import fabs
import Bio
from Bio.SeqRecord import SeqRecord
import re
from dataclasses import dataclass
from math import fabs



# Inheritance
@dataclass
class SearchSequnce:
    '''Finds sequence location from search query'''
    pass

@dataclass
class CRISPRCutting(SearchSequnce):
    '''Cuts DNA through CRISPR-cas9 double-stranded break.'''
    pass

@dataclass
class FindUpstreamAndDownStreamSequence(SearchSequnce):
    """Finds upstream and downstream sequences 
    from a sequence input in a search sequence."""
    pass

@dataclass
class CasEmbler(SearchSequnce): 
    """Simulate in vivo assembly and integration with the possibility"""
    pass


# Composition 
@dataclass
class SearchSequnce:
    '''Finds sequence location from search query'''
    sequence: Bio.SeqRecord
    chromosome: Bio.SeqRecord

    def __post_init__(self) -> None:
        self.start, self.end  , self.strand = find_sequence_location(self.sequence, self.chromosome)
    pass

@dataclass
class CRISPRCutting:
    '''Cuts DNA through CRISPR-cas9 double-stranded break.'''
    sequence: Bio.SeqRecord
    chromosome: Bio.SeqRecord
    def __post_init__(self) -> None:
        self.search_sequnce_obj = SearchSequnce(self.sequence,self.chromosome)
    pass

@dataclass
class FindUpstreamAndDownStreamSequence(SearchSequnce):
    """Finds upstream and downstream sequences 
    from a sequence input in a search sequence."""
    pass

@dataclass
class CasEmbler(SearchSequnce): 
    """Simulate in vivo assembly and integration with the possibility"""
    pass


### First iteration
@dataclass
class CRISPRCutting:
    """Cuts DNA through CRISPR-cas9 double-stranded break.
    
    Parameters
    ----------
    sequence : str
        This can be any search quury 
    chromosome : Bio.SeqRecord
        This could be any sequence to search in.
    """

    sequence: Bio.SeqRecord
    chromosome: Bio.SeqRecord

    def __post_init__(self) -> None:
        if len(self.sequence) != 20:
            raise ValueError(f"sgRNA must be 20 basepairs and not {len(self.sequence)}")

        # Check that we dont have multiple sgRNA sites
        self.find_all_occurences_of_a_sequence = find_all_occurences_of_a_sequence(
            self.sequence, self.chromosome
        )
        if self.find_all_occurences_of_a_sequence != 1:
            raise ValueError(
                f"sgRNA sites occuring in the chromosome should be 1 and not {len(self.find_all_occurences_of_a_sequence)}"
            )

        # Find searched sequence
        self.start_end_strand_index = find_sequence_location(
            self.sequence, self.chromosome
        )

        # Add attributes
        self.start_location = self.start_end_strand_index[0]
        self.end_location = self.start_end_strand_index[1]
        self.strand = self.start_end_strand_index[2]

        # CRISPR location
        self.crispr_base17_index = crispr_db_break_location(
            self.start_location, self.end_location, self.strand
        )
        self.crispr_seq = SeqRecord(
            self.chromosome.seq[self.start_location : self.end_location],
            name=f"sgRNA_{self.sequence.name}",
            id="",
        )

        # get upstream and sequences - get abosolut values and add annotations
        self.upstream_crispr_sequence = self.chromosome[
            : int(fabs(self.crispr_base17_index))
        ]
        self.downstream_crispr_sequence = self.chromosome[
            int(fabs(self.crispr_base17_index)) :
        ]

    def get_upstream_and_downstream_sequences(self):
        return self.upstream_crispr_sequence, self.downstream_crispr_sequence

    def add_annotation_to_upstream_and_downstream_sequences(self) -> None:
        # add annotations
        add_feature_annotation_to_seqrecord(
            self.upstream_crispr_sequence,
            label=f"{self.chromosome.name}_upstream_sequence_of_CRISPR_cut",
            strand=self.strand,
        )
        add_feature_annotation_to_seqrecord(
            self.downstream_crispr_sequence,
            label=f"{self.chromosome.name}_downstream_sequence_of_CRISPR_cut",
            strand=self.strand,
        )

    def display_cut(self):
        break_sequence = f"{str(self.upstream_crispr_sequence.seq)[-20:]}---x----{str(self.downstream_crispr_sequence.seq[:20])}"
        return break_sequence


@dataclass
class FindUpstreamAndDownStreamSequence:
    """Finds upstream and downstream sequences 
    from a sequence input in a search sequence.
    
    Parameters
    ----------
    sequence : str
        This can be any search quury 
    chromosome : Bio.SeqRecord
        This could be any sequence to search in.
    """

    sequence: Bio.SeqRecord
    chromosome: Bio.SeqRecord

    def __post_init__(self) -> None:
        # Find searched sequence
        self.start_end_strand_index = find_sequence_location(
            self.sequence, self.chromosome
        )

        # Add attributes
        self.start_location = self.start_end_strand_index[0]
        self.end_location = self.start_end_strand_index[1]
        self.strand = self.start_end_strand_index[2]

        # get upstream and sequences - get abosolut values and add annotations
        if self.strand > 0:
            self.upstream_sequence = self.chromosome[: int(fabs(self.start_location))]
            self.downstream_sequence = self.chromosome[int(fabs(self.end_location)) :]
        else:
            self.upstream_sequence = self.chromosome[: int(fabs(self.end_location))]
            self.downstream_sequence = self.chromosome[int(fabs(self.start_location)) :]

    def get_upstream_and_downstream_sequences(self):
        self.upstream_sequence.name, self.upstream_sequence.id = (
            f"Upstream_sequence_{self.chromosome.name}",
            f"Upstream_sequence_{self.chromosome.name}",
        )
        self.downstream_sequence.name, self.downstream_sequence.id = (
            f"Downstream_sequence_{self.chromosome.name}",
            f"Downstream_sequence_{self.chromosome.id}",
        )
        return self.upstream_sequence, self.downstream_sequence

    def add_annotation_to_upstream_and_downstream_sequences(self) -> None:
        # add annotations
        add_feature_annotation_to_seqrecord(
            self.upstream_sequence,
            label=f"{self.chromosome.name}_upstream_sequence",
            strand=self.strand,
        )
        add_feature_annotation_to_seqrecord(
            self.downstream_sequence,
            label=f"{self.chromosome.name}_downstream_sequence",
            strand=self.strand,
        )

    def find_occurrence_of_sequence_in_target_sequence(self):
        return find_all_occurences_of_a_sequence(self.sequence, self.chromosome)


def find_sequence_location(
    sequence: Bio.SeqRecord, sequence_to_search_in: Bio.SeqRecord
) -> tuple:
    """Finds start and end location of a mathced sequence.
    
    Parameters
    ----------
    sequence : str
    sequence_to_search_in : Bio.SeqRecord

    Returns
    -------
    (start_index,end_index) : tuple
    """
    strand = +1
    start_index = sequence_to_search_in.seq.find(sequence.seq)
    end_index = start_index + len(sequence)

    if start_index == -1:
        # search reverse_comp
        start_index = len(
            sequence_to_search_in
        ) - sequence_to_search_in.reverse_complement().seq.find(sequence.seq)
        end_index = start_index - len(sequence)
        strand = -1

        if start_index == -1:
            raise ValueError("ValueERROR - couldnt find a match")

    return (start_index, end_index, strand)


def add_feature_annotation_to_seqrecord(
    sequence: Bio.SeqRecord, label="", type_name="misc_feature", strand=0
) -> None:
    """ Adds feature, label and name to a Bio.Seqrecord sequence.
    Parameters
    ----------
    sequence : Bio.SeqRecord
    label : str (optional)
    type_name : str (default: "misc_feature")
    strand : int (default 0)
    
    Returns
    -------
    None
    """
    bio_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(sequence)), type=type_name, strand=strand
    )

    # label
    sequence.features.append(bio_feature)
    sequence.features[0].qualifiers["label"] = label
    sequence.features[0].qualifiers["name"] = sequence.name


def crispr_db_break_location(start_location, end_location, strand):
    """
    Determine the CRISPR cut location in the genome.

    Parameters
    ----------
    start_location : int
        Start position of the sgRNA sequence in the chromosome.
    end_location : int
        End position of the sgRNA sequence in the chromosome.
    strand : int
        Strand of the sgRNA sequence in the chromosome, +1 for positive strand, -1 for negative strand.

    Returns
    -------
    crispr_db_break : int
        CRISPR cut location in the genome.
    """
    if strand == +1:
        crispr_db_break = start_location + 17
    if strand == -1:
        crispr_db_break = start_location - 3

    return crispr_db_break


def seq_to_annotation(
    seq_record_from: Bio.SeqRecord, seq_record_onto: Bio.SeqRecord, type_name: str
):
    """Anotate an Bio.SeqRecord object from another amplicon object.

    Parameters
    ----------
    seqrec_from: Bio.SeqRecord
        annotation sequence that will be extracted

    seqrec_onto: Bio.SeqRecord

    type_name: str
        name of the sequence that will be extracted

    Returns
    -------
    None 
    """
    match_index = find_sequence_location(seq_record_from, seq_record_onto)
    start_location, end_location, strand = (
        match_index[0],
        match_index[1],
        match_index[2],
    )

    feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(start_location, end_location),
        type=type_name,
        strand=strand,
    )

    feature.qualifiers["label"] = seq_record_from.id
    seq_record_onto.features.append(feature)


def find_all_occurences_of_a_sequence(
    sequence: Bio.SeqRecord, sequence_to_search_in: Bio.SeqRecord
) -> tuple:
    """
    Searches for all occurrences of a given sequence in a given string.

    Parameters
    ----------
    sequence : Bio.SeqRecord
        Sequence to search for.
    sequence_to_search_in : Bio.SeqRecord
        Sequence to search in.

    Returns
    -------
    tuple
        Number of occurrences of `sequence` in `sequence_to_search_in`.

    """
    finder = re.finditer(
        str(sequence.seq.upper()), str(sequence_to_search_in.seq.upper())
    )
    matches_watson = [(match.start(), match.end()) for match in finder]

    if len(matches_watson) < 2:
        finder = re.finditer(
            str(sequence.seq.upper()),
            str(sequence_to_search_in.seq.reverse_complement().upper()),
        )
        matches_crick = [(match.start(), match.end()) for match in finder]

    return len(matches_watson) + len(matches_crick)
