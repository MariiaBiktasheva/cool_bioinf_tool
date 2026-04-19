from Bio import SeqIO
import os
from Bio.SeqUtils import gc_fraction

def filter_fastq (
    file_name: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
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
    
    #create output file name
    output_fastq = 'output_' + str(file_name.split('/')[1])
    
    #create folder for output data is it is necessary
    output_directory = "filtered"
    os.makedirs(output_directory, exist_ok=True)
    
   # open output file    
    with open(os.path.join(output_directory, output_fastq), "w") as output_handle:
        # process each record in input file
        for record in SeqIO.parse(file_name, "fastq"):
            
            # calculate gc content in percents
            gc = gc_fraction(record.seq)*100
            #calculate mean quality
            q = sum(record.letter_annotations['phred_quality']) / len(record.letter_annotations['phred_quality'])
            # calculate lengtha of the record
            l = len(record.seq)
            
            # check tresholds and bounds
            if  (l < length_bounds[1] and l > length_bounds[0] and 
                gc < gc_bounds[1] and gc > gc_bounds[0] and
                q > quality_threshold):
                
                # add  record to our output file
                SeqIO.write(record, output_handle, "fastq")
