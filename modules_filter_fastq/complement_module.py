# Add complement
def complement(sequence):
    """
    Returns a complementary sequence: attention!!!
    Returns only DNA
    """
    new_seq = ""
    Complement = {
        "A": "T",
        "a": "t",
        "T": "A",
        "t": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
    }

    # Select a complementary one for each nucleotide
    for nucl in sequence:
        new_seq=''.join([Complement[nucl] for nucl in sequence])
    return new_seq
