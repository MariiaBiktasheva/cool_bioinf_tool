from abc import ABC, abstractmethod

############################ Task. 1 Class  for sequences ############################


class BiologicalSequence(ABC):
    """This class gives interface to biological seqience data in general.
    It gives the possibilty to get length of sequence, to slice it, to print it in beautiful way

    Attributes:
    self.seq str -> sequence
    self.length  int -> length of sequence

    """

    def __init__(self, seq):

        self.seq = seq
        self.length = len(seq)

        # perform check of sequence just after the its initialization
        if not self.check_alphabet():
            raise ValueError("Invalid alphabet in sequence")

    def __getitem__(self, slc):
        return self.seq[slc]

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f"Sequence: {self.seq}\nSequence length: {self.length}"

    # abstract method of class for all bioseqs
    @abstractmethod
    def check_alphabet(self):
        """Check if seq is correct"""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for all nucleic acid sequences

    It performs the manipulations with seqs

    Methods:
    1. reverse
    2. coplement
    3. reverse-complement

    """

    @abstractmethod
    def complement(self):
        pass

    def reverse(self):
        return type(self)(self.seq[::-1])

    def reverse_complement(self):
        return self.reverse().complement()


class AminoAcidSequence(BiologicalSequence):
    """
    Class for all amino acid sequences

    It performs the manipulation with amino acid seqs

    Methods:
    1. get_aa_abundance: show abundances of aminoacid residuals in aa sequence
    2. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with one-letter aminoacid code
    """

    possible_letters = (
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    )

    def check_alphabet(self):
        return set(self.seq) <= set(self.possible_letters)

    def get_aa_abundance(self):
        aas = {}
        for aa in self.seq:
            if aa in aas:
                aas[aa] += 1
            else:
                aas[aa] = 1
        print(f"Aminoacid abundances:\n {aas}")


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences

    It performs the manipulation with DNA seqs

    Methods:
    1. trenscribe: convert dna to rna,  returns RNASequence object
    2. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with all dna- nucleotides
    """

    possible_letters = ("A", "T", "G", "C")

    def check_alphabet(self):
        return set(self.seq) <= set(self.possible_letters)

    def transcribe(self):
        rna = self.seq.replace("T", "U")
        return RNASequence(rna)

    def complement(self):
        return type(self)(self.seq.translate(str.maketrans("ATGC", "TACG")))


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences

    It performs the manipulation with DNA seqs

    Method:
    1. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with all rna- nucleotides
    """

    possible_letters = ("A", "U", "G", "C")

    def check_alphabet(self):
        return set(self.seq) <= set(self.possible_letters)

    def complement(self):
        return type(self)(self.seq.translate(str.maketrans("AGCU", "UCGA")))


############################ Task. 2 Class  for sequences ############################

from Bio import SeqIO
import os
from Bio.SeqUtils import gc_fraction


def filter_fastq(
    file_name: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filter sequences in fastq file by quality, length, GC content
    and write correct sequences into output file

    Args:
        file_name: fastq file with sequences
        quality_threshold: treshold for filtration by quality
        length_bounds: bounds of length to filter
        gc_bounds: bounds of GC content to filter

    Return:
        None (all correct sequence are written into output_file_name.fastq)
    """

    # process bounds
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    # create output file name
    output_fastq = "output_" + str(file_name.split("/")[1])

    # create folder for output data is it is necessary
    output_directory = "filtered"
    os.makedirs(output_directory, exist_ok=True)

    # open output file
    with open(os.path.join(output_directory, output_fastq), "w") as output_handle:
        # process each record in input file
        for record in SeqIO.parse(file_name, "fastq"):

            # calculate gc content in percents
            gc = gc_fraction(record.seq) * 100
            # calculate mean quality
            q = sum(record.letter_annotations["phred_quality"]) / len(
                record.letter_annotations["phred_quality"]
            )
            # calculate lengtha of the record
            l = len(record.seq)

            # check tresholds and bounds
            if (
                l < length_bounds[1]
                and l > length_bounds[0]
                and gc < gc_bounds[1]
                and gc > gc_bounds[0]
                and q > quality_threshold
            ):

                # add  record to our output file
                SeqIO.write(record, output_handle, "fastq")
