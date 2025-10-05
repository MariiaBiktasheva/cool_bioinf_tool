def gc_bounds_fun(seqs,
                  lower_bound: float = 0,
                  upper_bound: float = 100):
    """
    Filter sequences by GC content
    Args: 
        seq: dictionary with sequenses
        upper_bound: upper bound of GC content to filter
        lower_bound: lower bound of GC content to filter
    Return:
        Returns new dictionary new_seqs filtred by GC content 
    """
    gc_new_seqs={}
    for name in seqs:
         sequence = seqs[name][0]
         counts = 0
         for nucleotide in sequence:
             if nucleotide in 'GgCc':
                counts+=1
         gc_content_real = counts/len(sequence)*100
         if  lower_bound <= gc_content_real <= upper_bound: 
            gc_new_seqs[name] = seqs[name]
    return gc_new_seqs
