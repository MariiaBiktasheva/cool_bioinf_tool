# cool_bioinf_tool
The tool implements two functions.  
The first is  *run_dna_rna_tools* working with nucleotide sequences and the second is *filter_fastq* to filter reads according to specified parameters.
# USAGE
## run_dna_rna_tools
`run_dna_rna_tools` it only works with DNA and RNA, otherwise it gives an error. It takes sequences as str, and the last argument is the name of the function

__Functions__:
1. `is_nucleic_acid` performs a check, if not DNA and RNA
2. `transcribe` returns transcribed sequence
3. `reverse` returns reverse sequence
4. `complement` returns complement sequence
5. `reverse_complement` returns reverse&complement sequence
   __Return__: string with new sequences

## filter_fastq
`filter_fastq` filters sequences by quality, length, GC content

  __Args__:
  
        seqs: dictionary with sequences  
        quality_threshold: threshold for filtration by quality, default=0   
        length_bounds: bounds of length to filter, default (0, 2**32)  
        gc_bounds: bounds of GC content to filter, default (0,100) 
        
  __Return__: 
  
        Returns new dictionary filt_seqs filtered by quality, length, GC content
# Installation
Download as ordinary python files
## someday it will be real cool bioinf tool

