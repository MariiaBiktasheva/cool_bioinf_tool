
from Bio import SeqIO
import os
from Bio.SeqUtils import gc_fraction
import argparse
from loguru import logger

logger.remove(0)    # Удалить Handler по-умолчанию

logs_format = '<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>'

logger.add(       # Добавляем Handler для записи всех логов в файл all_logs.log
        "all_logs.log",
        level=0,
        format=logs_format,
        colorize=False,
        rotation="1 week"
    )
logger.add(       # Добавляем Handler для записи ошибок в файл errors_logs.log
        "errors_logs.log",
        level=29,
        format=logs_format,
        colorize=False,
        rotation="500 B"
    )


def filter_fastq (
    file_name: str,
    gc_bounds: float | list = [0, 100],
    length_bounds: int | list = [0, 2**32],
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
    if file_name is None:
        logger.error("No input file")
        raise TypeError

    elif not os.path.exists(file_name):
        logger.error("Wrong input path")
        raise ValueError

    # process bounds
    if gc_bounds is None:
        gc_bounds = (0, 100)
    elif isinstance(gc_bounds, list):
        if len(gc_bounds) == 1:
            gc_bounds = (0, gc_bounds[0])
        if len(gc_bounds) > 2:
            logger.error("More than 2 GC bounds")
            raise ValueError

    if length_bounds is None:
        length_bounds = (0, 2**32)

    elif isinstance(length_bounds, list):
        if len(length_bounds) == 1:
            length_bounds = (0, length_bounds[0])
        if len(length_bounds) > 2:
            logger.error("More than 2 length bounds : 2 first wil be used")
            raise ValueError

    if quality_threshold is None:
        quality_threshold = 0


    #create output file name
    output_fastq = 'output_' + str(file_name.split('/')[1])

    #create folder for output data is it is necessary
    output_directory = "filtered"
    os.makedirs(output_directory, exist_ok=True)

   # open output file    
    with open(os.path.join(output_directory, output_fastq), "w") as output_handle:
        # process each record in input file
        count=0
        count_total =0
        for record in SeqIO.parse(file_name, "fastq"):
            count_total +=1
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
                count+=1
                logger.info(f'{record.id} was successfully written to output file')

        logger.info(f'{count} sequensces out of {count_total} was successfully written to output file')

parser = argparse.ArgumentParser(
                        prog='Filter FASTQ files function',
                        description='This tool filters sequences in fastq file by quality, length, GC content and write correct sequences into output file ',
                        epilog='Text at the bottom of help')


parser.add_argument('-f', '--file_name', type=str, help='path to input fastq file (str)')
parser.add_argument('-gc', '--gc_bounds',type=float, nargs ='*', help='bounds of GC content to filter (float | list = [0, 100])')
parser.add_argument('-q', '--quality_threshold', type=float,  help='treshold for filtration by quality action (int = 0)')
parser.add_argument('-l', "--length_bounds", type=int, nargs ='*', help='bounds of length to filter (int | list = [0, 2**32])') 


if __name__ == '__main__':
    args = parser.parse_args()
    filter_fastq(file_name =args.file_name,
                 gc_bounds=args.gc_bounds,
                 length_bounds = args.length_bounds,
                 quality_threshold = args.quality_threshold)

