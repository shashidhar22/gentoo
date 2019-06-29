import os
import glob
import shutil
import warnings
import unittest
from gentoo import Index

class TestIndex(unittest.TestCase):

    def ignore_warnings(test_func):
        def do_test(self, *args, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                test_func(self, *args, **kwargs)
        return do_test

    def setUp(self):
        #Create test file and test folders
        dir_path = os.path.dirname(os.path.abspath(__file__))
        test_path = '{0}/tests'.format(dir_path)
        if not os.path.exists('tests'):
            os.mkdir(test_path)
        fastq_path = '{0}/Fastq'.format(test_path)
        dkc_path = '{0}/Dkc'.format(test_path)
        if not os.path.exists(dkc_path):
            os.mkdir(dkc_path)
        if not os.path.exists(dkc_path):
            os.mkdir(dkc_path)
        study_file = open('{0}/study_file.tsv'.format(test_path), 'w')
        sras = ['SRR6463548', 'SRR6463549']
        for sra in sras:
            study_file.write('{0}\t{0}\n'.format(sra))
        study_file.close()
        study_file = open('{0}/study_file2.tsv'.format(test_path), 'w')
        for sra in sras:
            study_file.write('{0}\t{1}/{0}_1.fastq.gz,{1}/{0}_2.fastq.gz\n'.format(sra, fastq_path))
        study_file.close()

    @ignore_warnings
    def test_groupFiles(self):
        file_list = ['test1.fas', 'test2.fa', 'test3.fasta',
                     'test4.fq', 'test5.fq.gz', 'test6.fastq', 'test7.fastq.gz',
                     'test8_r1.fq', 'test8_r2.fq', 'test9_R1.fq', 'test9_R2.fq',
                     'test10_r1.fq.gz', 'test10_r2.fq.gz', 'test11_R1.fq.gz', 'test11_R2.fq.gz',
                     'test12_l001_r1.fq', 'test12_l002_r1.fq', 'test12_l001_r2.fq', 'test12_l002_r2.fq',
                     'test13_l001_r1.fq.gz', 'test13_l002_r1.fq.gz', 'test13_l001_r2.fq.gz', 'test13_l002_r2.fq.gz',
                     'test14_l001_r1.fastq', 'test14_l002_r1.fastq', 'test14_l001_r2.fastq', 'test14_l002_r2.fastq',
                     'test15_l001_r1.fastq.gz', 'test15_l002_r1.fastq.gz', 'test15_l001_r2.fastq.gz', 'test15_l002_r2.fastq.gz',
                     'test16_l001_r1.fq', 'test16_l002_r1.fq', 'test16_l001_r2.fq', 'test16_l002_r2.fq',
                     'test16_l003_r1.fq', 'test16_l004_r1.fq', 'test16_l003_r2.fq', 'test16_l004_r2.fq',
                     'test17_l001_r1.fq.gz', 'test17_l002_r1.fq.gz', 'test17_l001_r2.fq.gz', 'test17_l002_r2.fq.gz',
                     'test17_l003_r1.fq.gz', 'test17_l004_r1.fq.gz', 'test17_l003_r2.fq.gz', 'test17_l004_r2.fq.gz',
                     'test18_l001_r1.fastq.gz', 'test18_l002_r1.fastq.gz', 'test18_l001_r2.fastq.gz', 'test18_l002_r2.fastq.gz',
                     'test18_l003_r1.fastq.gz', 'test18_l004_r1.fastq.gz', 'test18_l003_r2.fastq.gz', 'test18_l004_r2.fastq.gz',
                     'test19_l001_r1.fastq', 'test19_l002_r1.fastq', 'test19_l001_r2.fastq', 'test19_l002_r2.fastq',
                     'test19_l003_r1.fastq', 'test19_l004_r1.fastq', 'test19_l003_r2.fastq', 'test19_l004_r2.fastq',
                     'test20_L001_R1.fq', 'test20_L002_R1.fq', 'test20_L001_R2.fq', 'test20_L002_R2.fq',
                     'test21_L001_R1.fq.gz', 'test21_L002_R1.fq.gz', 'test21_L001_R2.fq.gz', 'test21_L002_R2.fq.gz',
                     'test22_L001_R1.fastq', 'test22_L002_R1.fastq', 'test22_L001_R2.fastq', 'test22_L002_R2.fastq',
                     'test23_L001_R1.fastq.gz', 'test23_L002_R1.fastq.gz', 'test23_L001_R2.fastq.gz', 'test23_L002_R2.fastq.gz',
                     'test24_L001_R1.fq', 'test24_L002_R1.fq', 'test24_L001_R2.fq', 'test24_L002_R2.fq',
                     'test24_L003_R1.fq', 'test24_L004_R1.fq', 'test24_L003_R2.fq', 'test24_L004_R2.fq',
                     'test25_L001_R1.fq.gz', 'test25_L002_R1.fq.gz', 'test25_L001_R2.fq.gz', 'test25_L002_R2.fq.gz',
                     'test25_L003_R1.fq.gz', 'test25_L004_R1.fq.gz', 'test25_L003_R2.fq.gz', 'test25_L004_R2.fq.gz',
                     'test26_L001_R1.fastq.gz', 'test26_L002_R1.fastq.gz', 'test26_L001_R2.fastq.gz', 'test26_L002_R2.fastq.gz',
                     'test26_L003_R1.fastq.gz', 'test26_L004_R1.fastq.gz', 'test26_L003_R2.fastq.gz', 'test26_L004_R2.fastq.gz',
                     'test27_L001_R1.fastq', 'test27_L002_R1.fastq', 'test27_L001_R2.fastq', 'test27_L002_R2.fastq',
                     'test27_L003_R1.fastq', 'test27_L004_R1.fastq', 'test27_L003_R2.fastq', 'test27_L004_R2.fastq',
                    ]
        indexer = Index('tests/', 'tests/', threads=2)
        study = indexer.groupFiles(file_list)
        true_study = {'test1': ['test1.fas'], 'test2': ['test2.fa'], 'test3': ['test3.fasta'],
                     'test4': ['test4.fq'], 'test5': ['test5.fq.gz'], 'test6': ['test6.fastq'], 'test7': ['test7.fastq.gz'],
                     'test8': ['test8_r1.fq', 'test8_r2.fq'], 'test9': ['test9_R1.fq', 'test9_R2.fq'],
                     'test10': ['test10_r1.fq.gz', 'test10_r2.fq.gz'], 'test11': ['test11_R1.fq.gz', 'test11_R2.fq.gz'],
                     'test12': ['test12_l001_r1.fq', 'test12_l002_r1.fq', 'test12_l001_r2.fq', 'test12_l002_r2.fq'],
                     'test13': ['test13_l001_r1.fq.gz', 'test13_l002_r1.fq.gz', 'test13_l001_r2.fq.gz', 'test13_l002_r2.fq.gz'],
                     'test14': ['test14_l001_r1.fastq', 'test14_l002_r1.fastq', 'test14_l001_r2.fastq', 'test14_l002_r2.fastq'],
                     'test15': ['test15_l001_r1.fastq.gz', 'test15_l002_r1.fastq.gz', 'test15_l001_r2.fastq.gz', 'test15_l002_r2.fastq.gz'],
                     'test16': ['test16_l001_r1.fq', 'test16_l002_r1.fq', 'test16_l001_r2.fq', 'test16_l002_r2.fq',
                                'test16_l003_r1.fq', 'test16_l004_r1.fq', 'test16_l003_r2.fq', 'test16_l004_r2.fq'],
                     'test17': ['test17_l001_r1.fq.gz', 'test17_l002_r1.fq.gz', 'test17_l001_r2.fq.gz', 'test17_l002_r2.fq.gz',
                                'test17_l003_r1.fq.gz', 'test17_l004_r1.fq.gz', 'test17_l003_r2.fq.gz', 'test17_l004_r2.fq.gz'],
                     'test18': ['test18_l001_r1.fastq.gz', 'test18_l002_r1.fastq.gz', 'test18_l001_r2.fastq.gz', 'test18_l002_r2.fastq.gz',
                                'test18_l003_r1.fastq.gz', 'test18_l004_r1.fastq.gz', 'test18_l003_r2.fastq.gz', 'test18_l004_r2.fastq.gz'],
                     'test19': ['test19_l001_r1.fastq', 'test19_l002_r1.fastq', 'test19_l001_r2.fastq', 'test19_l002_r2.fastq',
                                'test19_l003_r1.fastq', 'test19_l004_r1.fastq', 'test19_l003_r2.fastq', 'test19_l004_r2.fastq'],
                     'test20': ['test20_L001_R1.fq', 'test20_L002_R1.fq', 'test20_L001_R2.fq', 'test20_L002_R2.fq'],
                     'test21': ['test21_L001_R1.fq.gz', 'test21_L002_R1.fq.gz', 'test21_L001_R2.fq.gz', 'test21_L002_R2.fq.gz'],
                     'test22': ['test22_L001_R1.fastq', 'test22_L002_R1.fastq', 'test22_L001_R2.fastq', 'test22_L002_R2.fastq'],
                     'test23': ['test23_L001_R1.fastq.gz', 'test23_L002_R1.fastq.gz', 'test23_L001_R2.fastq.gz', 'test23_L002_R2.fastq.gz'],
                     'test24': ['test24_L001_R1.fq', 'test24_L002_R1.fq', 'test24_L001_R2.fq', 'test24_L002_R2.fq',
                                'test24_L003_R1.fq', 'test24_L004_R1.fq', 'test24_L003_R2.fq', 'test24_L004_R2.fq'],
                     'test25': ['test25_L001_R1.fq.gz', 'test25_L002_R1.fq.gz', 'test25_L001_R2.fq.gz', 'test25_L002_R2.fq.gz',
                                'test25_L003_R1.fq.gz', 'test25_L004_R1.fq.gz', 'test25_L003_R2.fq.gz', 'test25_L004_R2.fq.gz'],
                     'test26': ['test26_L001_R1.fastq.gz', 'test26_L002_R1.fastq.gz', 'test26_L001_R2.fastq.gz', 'test26_L002_R2.fastq.gz',
                                'test26_L003_R1.fastq.gz', 'test26_L004_R1.fastq.gz', 'test26_L003_R2.fastq.gz', 'test26_L004_R2.fastq.gz'],
                     'test27': ['test27_L001_R1.fastq', 'test27_L002_R1.fastq', 'test27_L001_R2.fastq', 'test27_L002_R2.fastq',
                                'test27_L003_R1.fastq', 'test27_L004_R1.fastq', 'test27_L003_R2.fastq', 'test27_L004_R2.fastq']
        }
        self.assertDictEqual(true_study, study, msg=[study, true_study])

    @ignore_warnings
    def test_sraDownload(self):
        true_sra = 'SRR6463548'
        indexer = Index('tests/', 'tests/')
        ret_code = indexer.sraDownload(true_sra)
        self.assertEqual(ret_code, 0)
        with self.assertRaises(SystemExit):
            indexer.sraDownload('ACC')

    @ignore_warnings
    def test_getStudy(self):
        indexer = Index('tests/study_file.tsv', 'tests/', threads=3)
        study = indexer.getStudy()
        true_study = dict()
        with open('tests/study_file.tsv') as study_file:
            for lines in study_file:
                lines = lines.strip().split('\t')[0]
                true_study[lines] = ['{0}/Fastq/{1}_1.fastq.gz'.format(indexer.out_path, lines), 
                                     '{0}/Fastq/{1}_2.fastq.gz'.format(indexer.out_path, lines)]
        self.assertDictEqual(true_study, study)

        indexer = Index('tests/study_file2.tsv', 'tests/', threads=1)
        study = indexer.getStudy()
        self.assertDictEqual(true_study, study)
        
    @ignore_warnings
    def test_prep(self):
        true_study = dict()
        out_path = os.path.abspath('tests/')
        with open('tests/study_file.tsv') as study_file:
            for lines in study_file:
                lines = lines.strip().split('\t')[0]
                true_study[lines] = ['{0}/Fastq/{1}_1.fastq.gz'.format(out_path, lines), 
                                     '{0}/Fastq/{1}_2.fastq.gz'.format(out_path, lines)]

        files = glob.glob('{0}/Fastq/*'.format(out_path))
        indexer = Index(files, 'tests/')
        study = indexer.prep()
        self.assertDictEqual(study, true_study)

        indexer = Index('tests/Fastq', 'tests/')
        study = indexer.prep()
        self.assertDictEqual(study, true_study)

        indexer = Index('tests/study_file.tsv', 'tests/', threads=3)
        study = indexer.prep()
        self.assertDictEqual(study, true_study)

        indexer = Index('tests/study_file2.tsv', 'tests/', threads=1)
        study = indexer.getStudy()
        self.assertDictEqual(true_study, study)

    def test_runKanalyze(self):
        sample = 'SRR6463548'
        out_dir = os.path.dirname(os.path.abspath(__file__))
        files = ['{0}/Fastq/{1}_1.fastq.gz'.format(out_dir, sample),
                '{0}/Fastq/{1}_2.fastq.gz'.format(out_dir, sample)]
        true_out = '{0}/tests/Dkc/{1}.dkc'.format(out_dir, sample)
        indexer = Index('tests/Fastq', 'tests/', threads=1, kmer=31)
        out_file, ret_code = indexer.runKanalyze((sample, files))
        self.assertEqual([true_out, 0], [out_file, ret_code])

    def test_createIndex(self):
        out_dir = os.path.dirname(os.path.abspath(__file__))
        true_out = [['{0}/tests/Dkc/{1}.dkc'.format(out_dir, 'SRR6463548'), 0],
                    ['{0}/tests/Dkc/{1}.dkc'.format(out_dir, 'SRR6463549'), 0]]
        samples = ['SRR6463548', 'SRR6463549']
        files = [['{0}/Fastq/{1}_1.fastq.gz'.format(out_dir, 'SRR6463548'), '{0}/Fastq/{1}_2.fastq.gz'.format(out_dir, 'SRR6463548')],
                 ['{0}/Fastq/{1}_1.fastq.gz'.format(out_dir, 'SRR6463549'), '{0}/Fastq/{1}_2.fastq.gz'.format(out_dir, 'SRR6463549')]]
        indexer = Index('tests/Fastq', 'tests/', temp_path='tests/')
        index_files = indexer.createIndex()
        self.assertEqual(true_out, index_files )




if __name__ == "__main__":
    unittest.main()