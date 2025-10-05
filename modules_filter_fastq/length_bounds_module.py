def length_bounds_fun(seqs, lower_bound: int = 0, upper_bound: int = 2**32):
    """
    filter sequences by length
    Args:
        seq: dictionary with sequenses
        upper_bound: upper bound of length to filter
        lower_bound: lower bound of length to filter
    Return:
        Returns new dictionary new_seqs filtred by length
    """
    len_new_seqs = {}
    for name in seqs:
        sequence = seqs[name][0]
        if lower_bound <= len(sequence) <= upper_bound:
            len_new_seqs[name] = seqs[name]
    return len_new_seqs
