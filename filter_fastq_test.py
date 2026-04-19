
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import unittest
import os
from  filter_fastq_oop import filter_fastq

# class for input checking
class TestInputErrors(unittest.TestCase):

    # check if empty file_name raise an error
    def test_input_empty(self):
        with self.assertRaises(TypeError):
            filter_fastq()

    # check if not existing path raise an error        
    def test_input_file_exists(self):
        with self.assertRaises(ValueError):
            filter_fastq(file_name = 'file/not/exists/')

    # check if 3 numbers in length bounds raise an error        
    def test_len_bounds(self):
        with self.assertRaises(ValueError):
            filter_fastq(file_name = 'data/example_fastq.fastq',  length_bounds = [10, 20, 30])

    # check if 3 number in gc bounds raise an error        
    def test_gc_bounds(self):
        with self.assertRaises(ValueError):
            filter_fastq(file_name = 'data/example_fastq.fastq',  gc_bounds = [10, 20, 30])


# class for output checking            
class TestOutPutFastq(unittest.TestCase):
    def setUp(self):
        self.inp = {'file_name' :'data/example_fastq.fastq', 
                    'length_bounds' : [10, 20],
                    'gc_bounds' : [30, 80],
                    'quality_threshold' : 37}

    def tearDown(self):
        try:
            os.remove('filtered/output_example_fastq.fastq')
        except FileNotFoundError:
            pass 

    # check if output file exists 
    def test_output_fastq_exists(self):
        filter_fastq(**self.inp)
        self.assertTrue(os.path.exists('filtered/output_example_fastq.fastq'))

    # check if output file is the same as correct test file      
    def test_output_fastq_file(self):
        filter_fastq(**self.inp)
        with open('filtered/output_example_fastq.fastq', 'r') as output:
            with open('data/test.fastq', 'r') as test:
                self.assertEqual(output.readlines(), test.readlines())

    # check if filtration is going well  
    def test_output_fastq_content(self):

        filter_fastq(**self.inp)

        for record in SeqIO.parse('filtered/output_example_fastq.fastq', "fastq"):

            # calculate gc content in percents
            gc = gc_fraction(record.seq)*100
            #calculate mean quality
            q = sum(record.letter_annotations['phred_quality']) / len(record.letter_annotations['phred_quality'])
            # calculate lengtha of the record
            l = len(record.seq)

            # check if length is ok 
            self.assertTrue(self.inp['length_bounds'][0] < l < self.inp['length_bounds'][1], 'Wrong length' )
             # check if gc is ok 
            self.assertTrue(self.inp['gc_bounds'][0] < gc < self.inp['gc_bounds'][1], 'Wrong gc' )
             # check if quality is ok 
            self.assertTrue( q > self.inp['quality_threshold'], 'Wrong quality' )

if __name__ == '__main__':
    unittest.main()



