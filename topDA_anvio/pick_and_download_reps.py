#!/usr/bin/env python3
import argparse, json, subprocess, sys, tempfile, shlex, os, pathlib, csv

# Rank helpers
REFSEQ_PREF = {"reference genome": 2, "representative genome": 1}
LEVEL_RANK  = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

def run(cmd):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr

def query_taxon(name, assembly_source):
    # Query NCBI Datasets for this exact species name (species rank only)
    # Then, if empty, try a broader text search
    base = f"datasets summary genome taxon {shlex.quote(name)} --assembly-version latest --as-json-lines --limit all --exclude-atypical --assembly-source {assembly_source} --tax-exact-match"
    rc, out, err = run(base)
    lines = [json.loads(l) for l in out.splitlines() if l.strip()]
    if not lines:
        # fallback: free-text search across species/assembly name/submitter
        fallback = f"datasets summary genome taxon {shlex.quote(name)} --assembly-version latest --as-json-lines --limit all --exclude-atypical --assembly-source {assembly_source} --search {shlex.quote(name)}"
        rc, out, err = run(fallback)
        lines = [json.loads(l) for l in out.splitlines() if l.strip()]
    return lines

def pick_one(cands):
    # Score each candidate:
    # 1) prefer refseq_category: reference > representative
    # 2) prefer assembly_level: Complete > Chromosome > Scaffold > Contig
    # 3) prefer RefSeq (GCF_) over GenBank (GCA_)
    # 4) prefer annotated assemblies
    # 5) newest release date
    best = None
    def score(rec):
        refcat = (rec.get("refseq_category") or "").lower()
        refcat_score = REFSEQ_PREF.get(refcat, 0)
        level = rec.get("assembly_level") or ""
        level_score = LEVEL_RANK.get(level, 0)
        acc = rec.get("accession") or rec.get("assembly", {}).get("accession") or ""
        is_refseq = 1 if acc.startswith("GCF_") else 0
        annotated = 1 if (rec.get("annotation_info") or {}).get("has_annotation") else 0
        date = rec.get("assembly", {}).get("release_date") or rec.get("release_date") or ""
        return (refcat_score, level_score, is_refseq, annotated, date)
    for r in cands:
        s = score(r)
        if best is None or s > best[0]:
            best = (s, r)
    return best[1] if best else None

def normalize(rec):
    # Flatten fields from the Genome Data Report
    asm = rec.get("assembly", {})
    org = rec.get("organism", {})
    ann = rec.get("annotation_info", {}) or {}
    return {
        "accession": rec.get("accession") or asm.get("accession"),
        "organism": org.get("organism_name") or rec.get("organism_name"),
        "refseq_category": (rec.get("refseq_category") or "").title(),
        "assembly_level": rec.get("assembly_level") or asm.get("assembly_level"),
        "source_db": ("RefSeq" if (rec.get("accession") or "").startswith("GCF_") else "GenBank" if (rec.get("accession") or "").startswith("GCA_") else ""),
        "release_date": asm.get("release_date") or rec.get("release_date") or "",
        "annotated": "yes" if ann.get("has_annotation") else "no",
        "assembly_name": asm.get("assembly_name") or rec.get("assembly_name") or "",
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species_file", help="Text file: one species name per line")
    ap.add_argument("--out", default="genomes_pkg", help="Output directory for the downloaded package")
    ap.add_argument("--include", default="genome", help="Comma-separated: genome,gff3,gbff,seq-report")
    args = ap.parse_args()

    with open(args.species_file) as f:
        names = [l.strip() for l in f if l.strip() and not l.strip().startswith("#")]

    selections, missing = [], []

    for name in names:
        # Try RefSeq first, then fall back to GenBank
        records = []
        for src in ("refseq", "all"):
            lines = query_taxon(name, src)
            # Filter out older versions if any 'is_latest' flag is present
            filtered = []
            for rec in lines:
                # Some schema variants put fields at top level, others under 'assembly'
                latest = rec.get("assembly", {}).get("is_latest") if isinstance(rec.get("assembly"), dict) else None
                # If the field is missing, we'll keep it (datasets already default to latest)
                if (latest is None) or (latest is True):
                    filtered.append(rec)
            records = filtered
            if records:
                break

        if not records:
            missing.append((name, "No assemblies found at NCBI"))
            continue

        chosen = pick_one(records)
        sel = normalize(chosen)
        sel["input_name"] = name
        selections.append(sel)

    # Write selection report
    pathlib.Path(args.out).mkdir(parents=True, exist_ok=True)
    rep_path = pathlib.Path("selection_report.tsv")
    with open(rep_path, "w", newline="") as tsv:
        w = csv.writer(tsv, delimiter="\t")
        w.writerow(["input_name","accession","organism","refseq_category","assembly_level","source_db","annotated","release_date","assembly_name"])
        for s in selections:
            w.writerow([s["input_name"], s["accession"], s["organism"], s["refseq_category"], s["assembly_level"], s["source_db"], s["annotated"], s["release_date"], s["assembly_name"]])
        for n,msg in missing:
            w.writerow([n, "", "", "", "", "", "", "", f"ERROR: {msg}"])

    # Collect accessions and download in one shot
    acc_file = pathlib.Path("chosen_accessions.txt")
    with open(acc_file, "w") as out:
        for s in selections:
            if s["accession"]:
                out.write(s["accession"] + "\n")

    if selections:
        zip_path = pathlib.Path("genomes.zip")
        cmd = f"datasets download genome accession --inputfile {shlex.quote(str(acc_file))} --include {shlex.quote(args.include)} --filename {shlex.quote(str(zip_path))}"
        rc, out, err = run(cmd)
        if rc != 0:
            print("ERROR: NCBI Datasets download failed:\n", err, file=sys.stderr)
            sys.exit(rc)
        # Unzip to output dir
        rc, out, err = run(f'unzip -q {shlex.quote(str(zip_path))} -d {shlex.quote(args.out)}')
        if rc != 0:
            print("ERROR: unzip failed:\n", err, file=sys.stderr)
            sys.exit(rc)

    # Print a concise status
    ok = len(selections)
    print(f"Selected {ok} assemblies. Details in selection_report.tsv; accessions in chosen_accessions.txt; files in {args.out}/")
    if missing:
        for n,msg in missing:
            print(f"[MISSING] {n}: {msg}", file=sys.stderr)

if __name__ == "__main__":
    main()

