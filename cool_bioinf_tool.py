# Import modules for filter_fastq
from modules_filter_fastq.gc_bounds_module import gc_bounds_fun
from modules_filter_fastq.length_bounds_module import length_bounds_fun
from modules_filter_fastq.quality_threshold_module import quality_threshold_fun

# Import modules for run_dna_rna_tools
from modules_filter_fastq.is_nucleic_acid_module import is_nucleic_acid
from modules_filter_fastq.transcribe_module import transcribe
from modules_filter_fastq.reverse_module import reverse
from modules_filter_fastq.complement_module import complement
from modules_filter_fastq.reverse_complement_module import reverse_complement


# Add filter_fastq
def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> dict[str, tuple[str, str]]:
    """
    Filter sequences by quality, length, GC content

    Args:
        seqs: dictionary with sequences
        quality_threshold: treshold for filtration by quality
        length_bounds: bounds of length to filter
        gc_bounds: bounds of GC content to filter

    Return:
        Returns new dictionary new_seqs filtered by quality, length, GC content
    """
    filt_seqs = {}
    filt_seqs = gc_bounds_fun(seqs, *gc_bounds)
    filt_seqs = length_bounds_fun(filt_seqs, *length_bounds)
    filt_seqs = quality_threshold_fun(filt_seqs, quality_threshold)
    return filt_seqs


# Add run_dna_rna_tools
def run_dna_rna_tools(*args) -> str:
    """
    Main function: it only works with DNA and RNA, otherwise it gives an error
    """
    # Setting the list for results
    results = []

    # Splitting arguments into sequences and name of function
    sequences = args[:-1]
    function = args[-1]

    # Creating a dictionary with names of functions
    tools = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    # Writing down function name
    function_name = tools[function]

    # Performing function
    for seq in sequences:
        if not is_nucleic_acid(seq):  # Performing a check, if not DNA and RNA
            print("Error:", seq, " is not DNA and RNA:((")
        else:
            result = function_name(seq)
            results.append(result)
    return " ".join(results)
