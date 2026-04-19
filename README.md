# it is HW21
# cool_bioinf_tool

The tool implements one functions.  

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
    
## Comand line usage example

`python filter_fastq_oop.py -f data/example_fastq.fastq -l 10 20`

## Python usage  example

`filter_fastq(file_name = 'data/example_fastq.fastq', 
                    length_bounds = [10, 20],
                    gc_bounds = [30, 80],
                    quality_threshold = 37)`



