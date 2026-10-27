#!/usr/bin/env python3
"""Read-only validation of authoritative per_genome.csv and per_inversion.csv.
Uses Python standard library only; does not repair or rewrite scientific data.
Usage: python scripts/audit_final_tables.py --genome path --inversion path --json report.json
"""
import argparse,csv,json,collections,pathlib,sys,hashlib
def audit(pg_path,pi_path):
  report={"inputs":{},"counts":{},"checks":{},"warnings":[]}
  for label,path in [("per_genome",pg_path),("per_inversion",pi_path)]:
    p=pathlib.Path(path)
    if not p.is_file(): raise FileNotFoundError(str(p))
    report["inputs"][label]={"file":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"bytes":p.stat().st_size}
  with open(pg_path,newline="",encoding="utf-8-sig") as f:
    pg=list(csv.DictReader(f))
  by_sample=collections.Counter(x["sample_id"] for x in pg)
  gk=collections.Counter(x["kingdom"] for x in pg)
  expected_k={"Metazoa":982,"Viridiplantae":181,"Fungi":11}
  report["counts"]["genomes"]=len(pg)
  report["counts"]["genomes_by_kingdom"]=dict(gk)
  report["checks"]["genomes_1174"]=len(pg)==1174
  report["checks"]["genomes_by_kingdom"]=dict(gk)==expected_k
  report["checks"]["unique_sample_ids"]=len(pg)==len(by_sample)
  reported={x["sample_id"]:int(x["inv_count_total"]) for x in pg}
  inv_ids=set(); duplicate_inv_ids=0; pi_counts=collections.Counter();ik=collections.Counter();total=0;below_200=0
  with open(pi_path,newline="",encoding="utf-8-sig") as f:
    reader=csv.DictReader(f)
    required={"inv_id","sample_id","kingdom","inv_len_bp"}
    if not required.issubset(reader.fieldnames or []): raise ValueError("missing inversion columns: "+str(required-set(reader.fieldnames or [])))
    for x in reader:
      total+=1
      sid=x["sample_id"]; pi_counts[sid]+=1; ik[x["kingdom"]]+=1
      if x["inv_id"] in inv_ids: duplicate_inv_ids+=1
      else: inv_ids.add(x["inv_id"])
      if float(x["inv_len_bp"])<200: below_200+=1
  report["counts"]["inversions"]=total
  report["counts"]["inversions_by_kingdom"]=dict(ik)
  report["counts"]["inversions_below_200bp"]=below_200
  report["counts"]["duplicate_inv_ids"]=duplicate_inv_ids
  report["checks"]["inversions_130159"]=total==130159
  report["checks"]["inversions_by_kingdom"]=dict(ik)=={"Metazoa":119022,"Viridiplantae":10976,"Fungi":161}
  report["checks"]["unique_inv_ids"]=duplicate_inv_ids==0
  report["checks"]["matched_sample_set"]=set(pi_counts)==set(by_sample)
  report["checks"]["genome_inv_counts_match"]=all(reported[s]==pi_counts[s] for s in reported)
  # A descriptive warning, not a mutation or invalidating assumption
  if below_200: report["warnings"].append(f"{below_200} inversion records are <200 bp; reconcile manuscript minimum-size threshold.")
  report["checks"]["all_required_checks"]=all(report["checks"].values())
  return report
def main():
  p=argparse.ArgumentParser()
  p.add_argument("--genome",default="5figure_Rplot/per_genome.csv")
  p.add_argument("--inversion",default="5figure_Rplot/per_inversion.csv")
  p.add_argument("--json",default="")
  a=p.parse_args(); result=audit(a.genome,a.inversion)
  content=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
  print(content)
  if a.json: pathlib.Path(a.json).write_text(content,encoding="utf-8")
  return 0 if result["checks"]["all_required_checks"] else 1
if __name__=="__main__": sys.exit(main())
