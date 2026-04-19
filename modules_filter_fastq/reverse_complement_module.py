# Add reverse_complement
def reverse_complement(sequence):
    """
    Returns a reverse complementary sequence:
    attention!!! It only returns DNA
    """
    return complement(reverse(sequence))
