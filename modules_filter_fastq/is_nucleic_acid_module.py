# Add is_nucleic_acid
def is_nucleic_acid(sequence) -> bool:
    """
    Check is sequence nucleic acid or no
    """
    # Setting the set of nucleotides separately for RNA and DNA
    rna_nucleotides = set("AaUuGgCc")
    dna_nucleotides = set("AaTtGgCc")

    # Check the nucleotides in sequence 
    unique_chars = set(sequence)
    is_dna = unique_chars <= dna_nucleotides
    is_rna = unique_chars <= rna_nucleotides
    mix = not (
        any(i in "Uu" for i in unique_chars)
        and any(i in "Tt" for i in unique_chars)
    )
    result = (is_dna or is_rna) and mix
    return result
