#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import unittest
## package
from leylab_pipelines.DB import Convert


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# global variables
ENSEMBL_ID = 'ENST00000407559'

HGNC_ID = '9245'
ENTREZ_ID = '39'
UNIPROT_ID = 'Q9BWD1'
ACCESSION_ID = 'AC131209'

EMAIL = 'dummy@dummy.info'
FAKE_EMAIL = 'pyentrez.info'


# tests
class Test_Convert(unittest.TestCase):

    def setUp(self):
        self.Id = Convert.Conversion(EMAIL)

    def tearDown(self):
        pass

    def test_email(self):
        self.assertRaises(ValueError, Convert.Conversion, FAKE_EMAIL)

    def test_convert_ensembl_to_entrez(self):
        # real ID
        ID = self.Id.convert_ensembl_to_entrez(ENSEMBL_ID)
        self.assertListEqual(ID, [ENSEMBL_ID, '55112'])
        # fake ID
        ID = self.Id.convert_ensembl_to_entrez(ENSEMBL_ID + '_FAKEID')
        self.assertListEqual(ID, [ENSEMBL_ID + '_FAKEID', None])
        
    def test_convert_hgnc_to_entrez(self):
        # real ID
        ID = self.Id.convert_hgnc_to_entrez(HGNC_ID)
        self.assertListEqual(ID, [HGNC_ID, '8500'])
        # fake ID
        ID = self.Id.convert_hgnc_to_entrez(HGNC_ID + '_FAKEID')
        self.assertListEqual(ID, [HGNC_ID + '_FAKEID', None])        

    def test_convert_entrez_to_uniprot(self):
        # real ID
        ID = self.Id.convert_entrez_to_uniprot(ENTREZ_ID)
        self.assertListEqual(ID, [ENTREZ_ID, 'Q9BWD1'])
        # fake ID

    def test_convert_uniprot_to_entrez(self):
        # real ID
        ID = self.Id.convert_uniprot_to_entrez(UNIPROT_ID)
        self.assertListEqual(ID, [UNIPROT_ID, '39'])
        # fake ID
        

    def test_convert_accesion_to_taxid(self):
        # real ID
        ID = self.Id.convert_accession_to_taxid(ACCESSION_ID)
        self.assertListEqual(ID, [ACCESSION_ID, '9606'])
        # fake ID
        ID = self.Id.convert_accession_to_taxid(ACCESSION_ID + '_FAKEID')
        self.assertListEqual(ID, [ACCESSION_ID + '_FAKEID', None])
    
