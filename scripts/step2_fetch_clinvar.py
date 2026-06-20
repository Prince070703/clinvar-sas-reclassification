import time
import os
import requests
import mysql.connector

# ================================================================
# STEP 2 — step2_fetch_clinvar.py  (FULLY FIXED)
#
# FIXES based on actual NCBI JSON structure:
#   1. Classification: germline_classification.description
#      (NOT clinical_significance — that field doesn't exist)
#   2. Disease: germline_classification.trait_set[].trait_name
#      (NOT top-level trait_set)
#   3. ref/alt: parsed from canonical_spdi (NC_...:pos:ref:alt)
#      because variation_loc always has ref="" alt=""
#   4. pos: variation_loc[].start (STRING — must cast to int)
#   5. HGVS: variation_set[].variation_name  (reliable field)
# ================================================================

DB = dict(host=os.getenv("MYSQL_HOST", "localhost"),
          user=os.getenv("MYSQL_USER", ""),
          password=os.getenv("MYSQL_PASSWORD", ""),
          database=os.getenv("MYSQL_DATABASE", "cancer_db"))

db   = mysql.connector.connect(**DB)
cur  = db.cursor(dictionary=True)
cur2 = db.cursor()

# ================================================================
# CONFIG
# ================================================================

NCBI_KEY  = os.getenv("NCBI_API_KEY", "")    # set in environment for faster NCBI requests
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RATE      = 0.12 if NCBI_KEY else 0.40
BATCH_DB  = 300
FETCH_SZ  = 200

GENES = [
    "BRCA1", "BRCA2", "PALB2", "ATM", "TP53",
    "PTEN",  "CDH1",  "APC",   "MLH1","MSH2",
]

CLASSIFICATIONS = [
    "Pathogenic",
    "Likely_pathogenic",
    "Uncertain_significance",
    "Likely_benign",
    "Benign",
    "Conflicting_interpretations_of_pathogenicity",
]

# ================================================================
# HELPERS
# ================================================================

def ncbi_get(endpoint, params):
    if NCBI_KEY:
        params["api_key"] = NCBI_KEY
    for attempt in range(4):
        try:
            r = requests.get(f"{NCBI_BASE}/{endpoint}",
                             params=params, timeout=45)
            if r.status_code == 200:
                return r
            print(f"      ⚠ HTTP {r.status_code} attempt {attempt+1}")
            time.sleep(2)
        except Exception as e:
            print(f"      ⚠ Network error attempt {attempt+1}: {e}")
            time.sleep(3)
    return None


def search_ids(gene, clnsig):
    """Return ALL ClinVar UIDs for gene + classification (paginated)."""
    all_ids  = []
    retmax   = 500
    retstart = 0

    while True:
        r = ncbi_get("esearch.fcgi", {
            "db":       "clinvar",
            "term":     f"{gene}[gene] AND {clnsig}[clinical significance]",
            "retmax":   retmax,
            "retstart": retstart,
            "retmode":  "json",
        })
        if not r:
            break

        data  = r.json().get("esearchresult", {})
        ids   = data.get("idlist", [])
        total = int(data.get("count", 0))
        all_ids.extend(ids)

        if len(all_ids) >= total or not ids:
            break
        retstart += retmax
        time.sleep(RATE)

    return all_ids


def parse_spdi(spdi):
    """
    Parse canonical_spdi to extract ref and alt alleles.

    Format: NC_000017.11:43091886:TTCTCTTCT:T
              accession : pos : deleted_seq : inserted_seq

    Returns (ref, alt) as strings, or (".", ".") if unparsable.
    """
    if not spdi:
        return ".", "."
    try:
        parts = spdi.split(":")
        if len(parts) == 4:
            ref = parts[2].strip() or "."
            alt = parts[3].strip() or "."
            # Cap at 500 chars (very long indels)
            return ref[:500], alt[:500]
    except Exception:
        pass
    return ".", "."


def fetch_details(uid_batch):
    """
    Fetch variant details for a batch of UIDs.
    Returns list of cleaned variant dicts.
    """
    if not uid_batch:
        return []

    r = ncbi_get("esummary.fcgi", {
        "db":      "clinvar",
        "id":      ",".join(str(u) for u in uid_batch),
        "retmode": "json",
    })
    if not r:
        return []

    try:
        data = r.json()
    except Exception as e:
        print(f"      ⚠ JSON parse error: {e}")
        return []

    result = data.get("result", {})
    uids_  = result.get("uids", [])
    rows   = []

    for uid in uids_:
        try:
            rec = result.get(str(uid), {})
            if not rec:
                continue

            # ── 1. CHROMOSOME + POSITION ──────────────────
            # From variation_set[0].variation_loc where assembly_name == "GRCh38"
            chrom = ""
            pos   = 0
            spdi  = ""

            for vset in rec.get("variation_set", []):

                # Get canonical_spdi for ref/alt parsing
                spdi = vset.get("canonical_spdi", "") or ""

                for loc in vset.get("variation_loc", []):
                    asm = loc.get("assembly_name", "")
                    # Match exactly "GRCh38" (confirmed from debug output)
                    if asm == "GRCh38":
                        chrom = str(loc.get("chr", "")).strip()
                        pos   = loc.get("start", "") or loc.get("display_start", "")
                        try:
                            pos = int(pos)
                        except (ValueError, TypeError):
                            pos = 0
                        break
                if chrom and pos:
                    break

            if not chrom or not pos:
                continue

            # Normalize chromosome
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            # ── 2. REF + ALT from canonical_spdi ─────────
            # NCBI variation_loc always has ref="" alt="" — must use SPDI
            ref, alt = parse_spdi(spdi)

            # ── 3. CLASSIFICATION ─────────────────────────
            # Key field: germline_classification (NOT clinical_significance)
            gc     = rec.get("germline_classification", {})
            clnsig = str(gc.get("description", "") or "").strip()
            review = str(gc.get("review_status", "") or "").strip()

            if not clnsig or clnsig.lower() in (
                    "not provided", "not classified", "other", ""):
                continue

            # Normalize to underscore format
            clnsig = clnsig.lower().replace(" ", "_").replace("/", "_")
            clnsig = clnsig[:100]

            # ── 4. DISEASE ────────────────────────────────
            # Under germline_classification.trait_set (NOT top-level)
            trait_set = gc.get("trait_set", [])
            disease   = "; ".join(
                str(t.get("trait_name", "") or "")
                for t in trait_set
                if t.get("trait_name")
            ) or "NA"
            disease = disease[:500]

            # ── 5. HGVS ──────────────────────────────────
            # variation_set[0].variation_name is the most reliable
            clnhgvs = ""
            for vset in rec.get("variation_set", []):
                name = str(vset.get("variation_name", "") or "")
                if name:
                    clnhgvs = name[:500]
                    break

            rows.append({
                "vid":     str(uid),
                "chrom":   chrom,
                "pos":     pos,
                "ref":     ref,
                "alt":     alt,
                "clnhgvs": clnhgvs,
                "clnsig":  clnsig,
                "review":  review[:200],
                "disease": disease,
            })

        except Exception as e:
            print(f"      ⚠ Parse error uid={uid}: {e}")
            continue

    return rows


def db_insert(gene, v):
    try:
        cur2.execute("""
            INSERT INTO clinvar_hboc10
                (chrom, pos, ref, alt, variant_id, gene,
                 clnhgvs, clinical_significance,
                 review_status, disease)
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ON DUPLICATE KEY UPDATE
                clinical_significance = VALUES(clinical_significance),
                review_status         = VALUES(review_status),
                disease               = VALUES(disease)
        """, (v["chrom"], v["pos"], v["ref"], v["alt"],
              v["vid"], gene, v["clnhgvs"],
              v["clnsig"], v["review"], v["disease"]))
        return True
    except Exception as e:
        print(f"      ⚠ DB error: {e}")
        return False


# ================================================================
# MAIN
# ================================================================

print("=" * 62)
print("  DISSERTATION PIPELINE — STEP 2 / 6")
print("  Fetching ClinVar via NCBI E-utilities API")
print(f"  Genes  : {GENES}")
print(f"  API key: {'YES' if NCBI_KEY else 'NO — add key for 3x speed'}")
print("=" * 62)

grand = 0

for gene in GENES:

    print(f"\n{'─'*62}")
    print(f"  Gene: {gene}")
    print(f"{'─'*62}")
    gtotal = 0

    for clnsig in CLASSIFICATIONS:

        print(f"  Searching: {clnsig:<50}", end="")
        uids = search_ids(gene, clnsig)
        print(f"{len(uids):>7} IDs")

        if not uids:
            time.sleep(RATE)
            continue

        batch_num = 0
        for start in range(0, len(uids), FETCH_SZ):
            batch_uids = uids[start : start + FETCH_SZ]
            batch_num += 1

            variants = fetch_details(batch_uids)
            inserted = 0

            for v in variants:
                if db_insert(gene, v):
                    inserted += 1
                    gtotal   += 1
                    grand    += 1

            if grand % BATCH_DB == 0 and grand > 0:
                db.commit()

            print(f"      batch {batch_num:>3}: "
                  f"{len(batch_uids):>3} IDs → "
                  f"{len(variants):>3} parsed → "
                  f"{inserted:>3} inserted  [total: {grand}]")

            time.sleep(RATE)

        db.commit()

    print(f"\n  ✔ {gene} DONE — {gtotal} variants")

# ── Final summary ─────────────────────────────────────────────
db.commit()
print("\n" + "=" * 62)
print("  FINAL SUMMARY")
print("=" * 62)
cur.execute("""
    SELECT gene,
           COUNT(*)  total,
           SUM(clinical_significance LIKE '%pathogenic%'
               AND clinical_significance NOT LIKE '%likely%'
               AND clinical_significance NOT LIKE '%conflict%') p,
           SUM(clinical_significance LIKE '%likely_pathogenic%')  lp,
           SUM(clinical_significance LIKE '%uncertain%')          vus,
           SUM(clinical_significance LIKE '%likely_benign%')      lb,
           SUM(clinical_significance LIKE '%benign%'
               AND clinical_significance NOT LIKE '%likely%')     b,
           SUM(clinical_significance LIKE '%conflict%')           conf
    FROM clinvar_hboc10
    GROUP BY gene ORDER BY gene
""")
print(f"\n  {'Gene':<8}{'Total':>7}{'P':>5}{'LP':>5}"
      f"{'VUS':>6}{'LB':>5}{'B':>5}{'Conf':>6}")
print(f"  {'-'*52}")
for r in cur.fetchall():
    print(f"  {r['gene']:<8}{r['total']:>7}{r['p']:>5}"
          f"{r['lp']:>5}{r['vus']:>6}{r['lb']:>5}"
          f"{r['b']:>5}{r['conf']:>6}")

db.close()
print(f"\n✅ STEP 2 DONE — {grand} total variants")
print("   Run next: python step3_fetch_gnomad.py\n")
