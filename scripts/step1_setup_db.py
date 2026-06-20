import os
import mysql.connector

# ================================================================
# step1_setup_db.py  — UPDATED SCHEMA (SAS-ONLY VERSION)
#
# CHANGES FROM PREVIOUS:
#   REMOVED: grpmax_faf, grpmax_pop columns entirely
#   KEPT:    ac_sas, an_sas, af_sas
#   ADDED:   faf95_sas  (Whiffin et al. 2017 formula)
#
# Only SAS population data is stored — no other populations
# ================================================================

DB = dict(host=os.getenv("MYSQL_HOST", "localhost"),
          user=os.getenv("MYSQL_USER", ""),
          password=os.getenv("MYSQL_PASSWORD", ""),
          database=os.getenv("MYSQL_DATABASE", "cancer_db"))

db  = mysql.connector.connect(**DB)
cur = db.cursor()

print("=" * 60)
print("  STEP 1 — Setting up database (SAS-only schema)")
print("=" * 60)

for tbl in ["reclassification_results",
            "gnomad_hboc10",
            "clinvar_hboc10"]:
    cur.execute(f"DROP TABLE IF EXISTS {tbl}")
    print(f"   Dropped: {tbl}")
db.commit()

# ── Table 1: ClinVar ──────────────────────────────────────────
cur.execute("""
CREATE TABLE clinvar_hboc10 (
    id                    INT AUTO_INCREMENT PRIMARY KEY,
    gene                  VARCHAR(20),
    chrom                 VARCHAR(10),
    pos                   INT,
    ref                   TEXT,
    alt                   TEXT,
    variant_id            VARCHAR(50),
    clnhgvs               TEXT,
    clinical_significance TEXT,
    review_status         TEXT,
    disease               TEXT,
    UNIQUE KEY uq_variant_gene (variant_id, gene),
    INDEX idx_gene (gene),
    INDEX idx_pos  (chrom, pos)
)
""")
print("   Created: clinvar_hboc10")

# ── Table 2: gnomAD — SAS ONLY, NO GRPMAX ─────────────────────
cur.execute("""
CREATE TABLE gnomad_hboc10 (
    id              INT AUTO_INCREMENT PRIMARY KEY,
    gene            VARCHAR(20),
    chrom           VARCHAR(10),
    pos             INT,
    ref             TEXT,
    alt             TEXT,

    ac_sas          INT   DEFAULT 0,
    an_sas          INT   DEFAULT 0,
    af_sas          FLOAT DEFAULT 0.0,
    faf95_sas       FLOAT DEFAULT 0.0,

    found_in_gnomad TINYINT DEFAULT 0,
    an_sufficient   TINYINT DEFAULT 0,

    INDEX idx_gene  (gene),
    INDEX idx_pos   (chrom, pos),
    INDEX idx_faf95 (faf95_sas),
    INDEX idx_found (found_in_gnomad)
)
""")
print("   Created: gnomad_hboc10  (SAS only — no grpmax)")

# ── Table 3: Results ──────────────────────────────────────────
cur.execute("""
CREATE TABLE reclassification_results (
    id                      INT AUTO_INCREMENT PRIMARY KEY,
    gene                    VARCHAR(20),
    chrom                   VARCHAR(10),
    pos                     INT,
    ref                     TEXT,
    alt                     TEXT,
    clnhgvs                 TEXT,
    disease                 TEXT,
    review_status           TEXT,
    original_classification TEXT,

    ac_sas                  INT   DEFAULT 0,
    an_sas                  INT   DEFAULT 0,
    af_sas                  FLOAT DEFAULT 0.0,
    faf95_sas               FLOAT DEFAULT 0.0,

    ba1_threshold           FLOAT,
    bs1_threshold           FLOAT,
    vcep_source             TEXT,

    criteria_applied        VARCHAR(30),
    new_classification      VARCHAR(50),
    change_direction        TEXT,

    found_in_gnomad         TINYINT DEFAULT 0,
    an_sufficient           TINYINT DEFAULT 0,

    INDEX idx_gene     (gene),
    INDEX idx_criteria (criteria_applied),
    INDEX idx_change   (change_direction(50))
)
""")
print("   Created: reclassification_results")

db.commit()
db.close()

print("""
======================================================
SCHEMA:
  ac_sas      = SAS allele count from gnomAD
  an_sas      = SAS allele number from gnomAD
  af_sas      = ac_sas / an_sas  (plain frequency)
  faf95_sas   = poisson.ppf(0.95, ac_sas) / an_sas
                (Whiffin et al. 2017 — 95% CI upper bound)

  NO grpmax_faf
  NO grpmax_pop
  NO other populations

Run next: python step2_fetch_clinvar.py
======================================================
""")
