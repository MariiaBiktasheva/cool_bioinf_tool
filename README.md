# cool_bioinf_tool

## it is HW16
## cool_bioinf_tool

The tool implements two functions.

## First is actions with biological sequences

## Class BiologicalSequence(ABC):
    This class gives interface to biological seqience data in general.
    It gives the possibilty to get length of sequence, to slice it, to print it in beautiful way

    Attributes:
    self.seq str -> sequence
    self.length  int -> length of sequence

## Class NucleicAcidSequence(BiologicalSequence):
    
    Class for all nucleic acid sequences

    It performs the manipulations with seqs

    Methods:
    1. reverse
    2. coplement
    3. reverse-complement

## Class AminoAcidSequence(BiologicalSequence):
    
    Class for all amino acid sequences

    It performs the manipulation with amino acid seqs

    Methods:
    1. get_aa_abundance: show abundances of aminoacid residuals in aa sequence
    2. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with one-letter aminoacid code

### Class DNASequence(NucleicAcidSequence):
   
    Class for DNA sequences

    It performs the manipulation with DNA seqs

    Methods:
    1. trenscribe: convert dna to rna,  returns RNASequence object
    2. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with all dna- nucleotides

### Class RNASequence(NucleicAcidSequence):
 
    Class for RNA sequences

    It performs the manipulation with RNA seqs

    Method:
    1. check_alphabet: check if seqs content only correct letters

    Atributes:
    1. possible_letters: tuple with all rna- nucleotides

## See usage examples in exaples.ipynb

## Second function is filtering fastq files


## filter_fastq
`filter_fastq` filters sequences by quality, length, GC content

  
    Filter sequences in fastq file by quality, length, GC content 
    and write correct sequences into output file 

    Args:
        file_name: fastq file with sequences
        quality_threshold: treshold for filtration by quality
        length_bounds: bounds of length to filter
        gc_bounds: bounds of GC content to filter

    Return:
        None (all correct sequence are written into output_file_name.fastq)
    

## Python usage  example

`filter_fastq(file_name = 'data/example_fastq.fastq', 
                    length_bounds = [10, 20],
                    gc_bounds = [30, 80],
                    quality_threshold = 37)`



