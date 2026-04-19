def quality_threshold_fun(seqs, quality_threshold: int = 0):
    """
    filters sequences by quality
    Args:
        seqs: dictionary with sequenses
        quality_threshold: treshold for filtration
    Return:
        Returns new dictionary new_seqs filtred by quality
    """
    q_new_seqs = {}
    for name in seqs:
        characters = seqs[name][1]
        sum_quality = 0
        for char in characters:
            sum_quality += ord(char) - 33
        mean_quality = sum_quality / len(characters)
        if quality_threshold < mean_quality:
            q_new_seqs[name] = seqs[name]
    return q_new_seqs
