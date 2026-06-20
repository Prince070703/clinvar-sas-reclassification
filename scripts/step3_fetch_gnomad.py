import time
import os
import requests
import mysql.connector
from scipy.stats import poisson

# ================================================================
# step3_fetch_gnomad.py  — SAS FAF95 ONLY (FIXED VERSION)
#
# FIXES FROM PREVIOUS VERSION:
#   1. REMOVED position fallback (by_pos)
#      Old: if no exact match → grab first gnomAD at same pos
#      New: if no exact match → mark as ABSENT (PM2_Supporting)
#      Why: assigning another variant's frequency is scientifically
#           wrong. Absence = PM2_Supporting which is correct ACMG.
#
#   2. EXACT MATCH ONLY: (chrom, pos, ref, alt) must all match
#      This prevents high-frequency benign variants being
#      incorrectly assigned to unrelated pathogenic entries.
#
#   3. gnomAD uses 1-based positions, NCBI uses 0-based
#      We try both pos and pos+1 for matching
#
# FAF95 formula — Whiffin et al. 2017:
#   FAF95 = poisson.ppf(0.95, AC_sas) / AN_sas
#   95% Poisson CI upper bound — more conservative than AC/AN
# ================================================================

DB   = dict(host=os.getenv("MYSQL_HOST", "localhost"),
            user=os.getenv("MYSQL_USER", ""),
            password=os.getenv("MYSQL_PASSWORD", ""),
            database=os.getenv("MYSQL_DATABASE", "cancer_db"))
db   = mysql.connector.connect(**DB)
cur  = db.cursor(dictionary=True)
cur2 = db.cursor()

URL        = "https://gnomad.broadinstitute.org/api"
MIN_AN_SAS = 2000

GENES = ["BRCA1","BRCA2","PALB2","ATM","TP53",
         "PTEN","CDH1","APC","MLH1","MSH2"]

QUERY = """
query ($gene: String!, $dataset: DatasetId!) {
  gene(gene_symbol: $gene, reference_genome: GRCh38) {
    variants(dataset: $dataset) {
      pos
      ref
      alt
      exome  { populations { id ac an } }
      genome { populations { id ac an } }
    }
  }
}
"""

# ================================================================
# FAF95 — SAS ONLY
# ================================================================

def calc_faf95(ac, an):
    """
    FAF95 = poisson.ppf(0.95, AC_sas) / AN_sas
    Whiffin et al. Genet Med 2017

    AC=0 or AN=0  → 0.0
    AC >= AN      → 1.0 (capped)
    Result > 1.0  → 1.0 (capped — edge case when AC is very
                    close to AN, Poisson upper CI slightly
                    exceeds AN giving FAF95 just above 1.0)
    """
    if not ac or not an or an == 0:
        return 0.0
    if ac >= an:
        return 1.0
    faf95 = float(poisson.ppf(0.95, ac)) / an
    return round(min(faf95, 1.0), 10)   # always cap at 1.0


def get_sas(populations):
    """
    Extract SAS AC and AN from gnomAD populations list.
    Returns (ac, an) for South Asian population only.
    Returns (0, 0) if SAS not found.
    """
    for p in populations:
        if str(p.get("id","")).lower() == "sas":
            ac = int(p.get("ac") or 0)
            an = int(p.get("an") or 0)
            return ac, an
    return 0, 0


def get_combined_sas(variant):
    """
    Combine SAS AC/AN across gnomAD exome and genome calls for the
    same variant, matching the dissertation wording: exomes + genomes.
    """
    total_ac = 0
    total_an = 0
    for source in ("exome", "genome"):
        data = variant.get(source) or {}
        ac, an = get_sas(data.get("populations", []))
        total_ac += ac
        total_an += an
    return total_ac, total_an


# ================================================================
# FETCH FROM gnomAD API
# ================================================================

def fetch_gene(gene):
    """Fetch all variants for a gene from gnomAD v4.1"""
    print(f"   Querying gnomAD: {gene}...", end=" ", flush=True)
    try:
        r = requests.post(
            URL,
            json={"query": QUERY,
                  "variables": {"gene": gene,
                                "dataset": "gnomad_r4"}},
            headers={"Content-Type": "application/json"},
            timeout=300
        )
        if r.status_code != 200:
            print(f"HTTP {r.status_code}")
            return []
        data = r.json()
        if data.get("errors"):
            print(f"Error: {data['errors'][0].get('message','')}")
            return []
        variants = (data.get("data",{})
                       .get("gene",{})
                       .get("variants",[]) or [])
        print(f"{len(variants)} variants fetched")
        return variants
    except Exception as e:
        print(f"FAILED — {e}")
        return []


# ================================================================
# INSERT HELPERS
# ================================================================

def insert_matched(gene, ch, pos, ref, alt, ac, an, af, f95, an_ok):
    """Insert a variant that was found in gnomAD SAS"""
    cur2.execute("""
        INSERT INTO gnomad_hboc10
            (gene, chrom, pos, ref, alt,
             ac_sas, an_sas, af_sas, faf95_sas,
             found_in_gnomad, an_sufficient)
        VALUES (%s,%s,%s,%s,%s,
                %s,%s,%s,%s, 1,%s)
    """, (gene, ch, pos, ref, alt,
          ac, an, af, f95, an_ok))


def insert_absent(gene, ch, pos, ref, alt):
    """Insert a variant NOT found in gnomAD — absent = PM2_Supporting"""
    cur2.execute("""
        INSERT INTO gnomad_hboc10
            (gene, chrom, pos, ref, alt,
             ac_sas, an_sas, af_sas, faf95_sas,
             found_in_gnomad, an_sufficient)
        VALUES (%s,%s,%s,%s,%s,
                0, 0, 0, 0, 0, 0)
    """, (gene, ch, pos, ref, alt))


# ================================================================
# MAIN
# ================================================================

print("=" * 62)
print("  STEP 3 — Fetch gnomAD SAS frequencies (FIXED)")
print("  Metric : FAF95_sas = poisson.ppf(0.95, AC_sas) / AN_sas")
print("  Match  : EXACT (chrom + pos + ref + alt) only")
print("  No pos fallback — unmatched = absent = PM2_Supporting")
print("  NO grpmax — SAS population ONLY")
print("=" * 62)

# Clear previous gnomAD data
cur2.execute("TRUNCATE TABLE gnomad_hboc10")
db.commit()
print("\n   Cleared gnomad_hboc10\n")

total_matched  = 0
total_absent   = 0
total_nospdi   = 0

for gene in GENES:
    print(f"\n{'─'*62}")
    print(f"  {gene}")
    print(f"{'─'*62}")

    # Get ClinVar variants for this gene
    cur.execute("""
        SELECT DISTINCT chrom, pos, ref, alt
        FROM clinvar_hboc10
        WHERE gene = %s
    """, (gene,))
    clinvars = cur.fetchall()
    print(f"   ClinVar variants : {len(clinvars):,}")

    # Fetch gnomAD variants for this gene
    gnomad = fetch_gene(gene)

    # ── Build EXACT lookup: (pos, ref, alt) → variant ─────────
    # gnomAD uses 1-based positions
    # NCBI/ClinVar uses 0-based start positions
    # We index BOTH pos and pos-1 to handle offset
    exact = {}
    for gv in gnomad:
        try:
            p   = int(gv["pos"])
            ref = str(gv["ref"] or "").upper()
            alt = str(gv["alt"] or "").upper()
            ac, an = get_combined_sas(gv)
            af     = round(ac/an, 10) if an > 0 else 0.0
            f95    = calc_faf95(ac, an)
            an_ok  = 1 if an >= MIN_AN_SAS else 0

            entry = (ac, an, af, f95, an_ok)

            # Index by both gnomAD pos and pos-1 (coordinate offset)
            exact[(p,   ref, alt)] = entry
            exact[(p-1, ref, alt)] = entry

        except Exception:
            continue

    matched = absent = nospdi = 0

    for cv in clinvars:
        pos = int(cv["pos"])
        ref = str(cv["ref"] or "").upper()
        alt = str(cv["alt"] or "").upper()
        ch  = cv["chrom"]

        # ── No SPDI structural variants ────────────────────────
        if ref in ("", ".") or alt in ("", "."):
            insert_absent(gene, ch, pos, ref, alt)
            nospdi += 1
            continue

        # ── EXACT MATCH ONLY ───────────────────────────────────
        entry = exact.get((pos, ref, alt))

        if entry:
            ac, an, af, f95, an_ok = entry
            insert_matched(gene, ch, pos, ref, alt,
                           ac, an, af, f95, an_ok)
            matched += 1
        else:
            # No exact match → ABSENT → PM2_Supporting
            insert_absent(gene, ch, pos, ref, alt)
            absent += 1

    db.commit()

    total_matched += matched
    total_absent  += absent
    total_nospdi  += nospdi

    # ── Show top 5 FAF95 values for this gene ──────────────────
    cur.execute("""
        SELECT ac_sas, an_sas, af_sas, faf95_sas
        FROM gnomad_hboc10
        WHERE gene = %s
          AND found_in_gnomad = 1
          AND faf95_sas > 0
          AND an_sas >= 2000
        ORDER BY faf95_sas DESC
        LIMIT 5
    """, (gene,))
    samples = cur.fetchall()

    print(f"   Matched (exact)  : {matched:,}")
    print(f"   Absent (PM2)     : {absent:,}")
    print(f"   No SPDI          : {nospdi:,}")

    if samples:
        print(f"   Top FAF95_sas values (AN >= 2000):")
        for s in samples:
            ac  = s['ac_sas']
            an  = s['an_sas']
            af  = float(s['af_sas'])
            f95 = float(s['faf95_sas'])
            print(f"     AC={ac:>6}  AN={an:>6}  "
                  f"AF={af:.6f}  FAF95={f95:.6f}")
    else:
        print(f"   No variants with AN >= 2000 found in SAS")

    time.sleep(3)

db.commit()

# ================================================================
# FINAL SUMMARY
# ================================================================

print(f"\n{'='*62}")
print(f"✅ STEP 3 COMPLETE")
print(f"\n   Matched (exact)  : {total_matched:,}")
print(f"   Absent (PM2)     : {total_absent:,}")
print(f"   No SPDI          : {total_nospdi:,}")

print(f"\n   Method: EXACT match on (pos, ref, alt)")
print(f"   Unmatched → absent → PM2_Supporting (correct ACMG)")
print(f"   SAS FAF95 = poisson.ppf(0.95, AC_sas) / AN_sas")
print(f"   No grpmax — SAS population only")

# Per gene summary
print(f"\n{'─'*62}")
print(f"   PER-GENE SUMMARY:")
cur.execute("""
    SELECT
        gene,
        COUNT(*)                    total,
        SUM(found_in_gnomad = 1)    found,
        SUM(found_in_gnomad = 0)    absent,
        SUM(an_sufficient = 1)      an_ok,
        SUM(faf95_sas > 0)          has_faf95,
        ROUND(MAX(faf95_sas), 6)    max_faf95,
        SUM(faf95_sas >= 0.001)     above_ba1_brca,
        SUM(faf95_sas >= 0.0001)    above_bs1_brca
    FROM gnomad_hboc10
    GROUP BY gene ORDER BY gene
""")

print(f"\n   {'Gene':<8}{'Total':>8}{'Found':>8}"
      f"{'Absent':>8}{'AN≥2k':>7}"
      f"{'MaxFAF95':>10}")
print(f"   {'─'*52}")

for r in cur.fetchall():
    print(f"   {r['gene']:<8}"
          f"{r['total']:>8,}"
          f"{r['found']:>8,}"
          f"{r['absent']:>8,}"
          f"{r['an_ok']:>7,}"
          f"{float(r['max_faf95'] or 0):>10.6f}")

# Sanity check
print(f"\n   SANITY CHECK — FAF95 distribution:")
cur.execute("""
    SELECT
        SUM(faf95_sas = 0)            AS zero,
        SUM(faf95_sas > 0
            AND faf95_sas < 0.0001)   AS very_rare,
        SUM(faf95_sas >= 0.0001
            AND faf95_sas < 0.001)    AS bs1_range,
        SUM(faf95_sas >= 0.001)       AS ba1_range,
        MAX(faf95_sas)                AS max_val,
        AVG(CASE WHEN faf95_sas > 0
                 THEN faf95_sas END)  AS avg_nonzero
    FROM gnomad_hboc10
    WHERE found_in_gnomad = 1
      AND an_sas >= 2000
""")
s = cur.fetchone()
if s:
    print(f"     FAF95 = 0          : {s['zero']:,}")
    print(f"     0 < FAF95 < 0.0001 : {s['very_rare']:,}  (very rare)")
    print(f"     0.0001-0.001       : {s['bs1_range']:,}  (BS1 range)")
    print(f"     >= 0.001           : {s['ba1_range']:,}  (BA1 range)")
    print(f"     Max FAF95          : {float(s['max_val'] or 0):.6f}")
    print(f"     Avg FAF95 (>0)     : {float(s['avg_nonzero'] or 0):.8f}")

db.close()
print(f"\n{'='*62}")
print("   Run next: python step4_reclassify.py\n")
