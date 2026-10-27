import csv, pathlib, tempfile, unittest, importlib.util
MODULE=pathlib.Path(__file__).resolve().parents[1]/"scripts"/"audit_final_tables.py"
spec=importlib.util.spec_from_file_location("audit_final_tables", MODULE)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
class TestAudit(unittest.TestCase):
 def test_fixture_expected_fail_on_tiny_counts(self):
  with tempfile.TemporaryDirectory() as td:
   pg=pathlib.Path(td)/"per_genome.csv"; pi=pathlib.Path(td)/"per_inversion.csv"
   pg.write_text("sample_id,kingdom,inv_count_total\ns1,Metazoa,1\n",encoding="utf-8")
   pi.write_text("inv_id,sample_id,kingdom,inv_len_bp\ni1,s1,Metazoa,220\n",encoding="utf-8")
   result=audit.audit(pg,pi)
   self.assertEqual(result["counts"]["inversions"],1)
   self.assertTrue(result["checks"]["unique_inv_ids"])
   self.assertTrue(result["checks"]["genome_inv_counts_match"])
   self.assertFalse(result["checks"]["genomes_1174"])
   self.assertFalse(result["checks"]["all_required_checks"])
 def test_duplicate_ids_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   pg=pathlib.Path(td)/"pg.csv";pi=pathlib.Path(td)/"pi.csv"
   pg.write_text("sample_id,kingdom,inv_count_total\ns1,Metazoa,2\n",encoding="utf-8")
   pi.write_text("inv_id,sample_id,kingdom,inv_len_bp\ni1,s1,Metazoa,199\ni1,s1,Metazoa,201\n",encoding="utf-8")
   r=audit.audit(pg,pi)
   self.assertEqual(r["counts"]["duplicate_inv_ids"],1)
   self.assertEqual(r["counts"]["inversions_below_200bp"],1)
   self.assertFalse(r["checks"]["unique_inv_ids"])
if __name__=="__main__": unittest.main()
