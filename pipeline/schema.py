"""SQLite schema. Spec Section 7, column for column.

Adapted from debris/pipeline/schema.py (see REUSE.md): the read-only connect and
the build-then-swap convention are kept, the tables are WHITTLE's.

The droplet gets this file and the gzipped structures, nothing else, and it does
no request-time computation (spec 7, last paragraph). Everything the front page
needs is precomputed into app/static/data/*.json.
"""
from __future__ import annotations

import sqlite3

from pipeline.common import DB_PATH, connect, log

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS targets (
    target              TEXT PRIMARY KEY,
    uniprot             TEXT,
    name                TEXT,
    construct           TEXT,
    range_start         INTEGER,
    range_end           INTEGER,
    chain_length        INTEGER,
    stoichiometry       INTEGER,
    his_enriched        INTEGER,
    ph_switch_intended  INTEGER,
    known_nanobodies    INTEGER,
    ecd_scoped          INTEGER,
    antigen_md5         TEXT,
    tier2_cost          TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS candidates (
    id                          TEXT PRIMARY KEY,
    target                      TEXT NOT NULL REFERENCES targets(target),
    sequence                    TEXT NOT NULL,
    seq_md5                     TEXT NOT NULL,
    length                      INTEGER,
    multi_target                INTEGER DEFAULT 0,
    numbering_method            TEXT,
    v_gene                      TEXT,
    cdr1                        TEXT,
    cdr2                        TEXT,
    cdr3                        TEXT,
    cdr1_len                    INTEGER,
    cdr2_len                    INTEGER,
    cdr3_len                    INTEGER,
    cdr3_cluster                INTEGER,
    hallmark_intact             INTEGER,   -- strict: all four camelid residues
    hallmark_recomputed         INTEGER,
    hallmark_n                  INTEGER,   -- 0-4, how many of the tetrad are camelid
    hallmark_detail             TEXT,
    his_cdr_total               INTEGER,
    his_cdr1                    INTEGER,
    his_cdr2                    INTEGER,
    his_cdr3                    INTEGER,
    net_charge_74               REAL,
    pi                          REAL,
    canonical_disulfide         INTEGER,
    extra_disulfide             INTEGER,
    n_liabilities               INTEGER DEFAULT 0,
    -- upstream columns, retained for comparison only (spec 5.1 and Appendix B item 1)
    upstream_tm_pred            REAL,
    upstream_liability_score    REAL,
    upstream_on_nll             REAL,
    upstream_score              REAL,
    upstream_s_tm               REAL,
    upstream_s_liab             REAL,
    upstream_s_nll              REAL,
    upstream_s_his              REAL,
    upstream_hallmark           INTEGER,
    upstream_cdr3_cluster       INTEGER,
    upstream_his_cdr_total      INTEGER,
    upstream_his_cdr1           INTEGER,
    upstream_his_cdr2           INTEGER,
    upstream_his_cdr3           INTEGER,
    upstream_target_cdr_identity   REAL,
    upstream_global_cdr_identity   REAL,
    upstream_global_nearest        TEXT,
    upstream_global_n_known        INTEGER,
    upstream_cdr_novelty        REAL,   -- empty in all five source files (spec 5.1 item 3)
    upstream_cdr_identity       REAL    -- ditto
);
CREATE INDEX IF NOT EXISTS idx_candidates_target ON candidates(target);
CREATE INDEX IF NOT EXISTS idx_candidates_md5 ON candidates(seq_md5);
CREATE INDEX IF NOT EXISTS idx_candidates_cluster ON candidates(target, cdr3_cluster);

CREATE TABLE IF NOT EXISTS liabilities (
    candidate_id   TEXT NOT NULL REFERENCES candidates(id),
    kind           TEXT NOT NULL,
    region         TEXT,
    position       INTEGER,
    detail         TEXT
);
CREATE INDEX IF NOT EXISTS idx_liabilities_cand ON liabilities(candidate_id);

CREATE TABLE IF NOT EXISTS structures (
    candidate_id         TEXT PRIMARY KEY REFERENCES candidates(id),
    path                 TEXT,
    refined              INTEGER,
    mean_error           REAL,
    cdr3_error           REAL,
    psh                  REAL,      -- TNP patch surface hydrophobicity
    ppc                  REAL,      -- patch positive charge
    pnc                  REAL,      -- patch negative charge
    sfvcsp_na            REAL,      -- SFvCSP: not applicable to a single-domain VHH
    tnp_flags            TEXT,
    tnp_version          TEXT,
    tnp_percentile       REAL,
    swi                  REAL,
    protein_sol          REAL,
    aggrescan3d          REAL,
    hydrophobic_patch_area REAL,
    charged_patch_area   REAL,
    sasa_total           REAL,
    nativeness_ppl       REAL,
    oasis_humanness      REAL,
    nearest_pdb          TEXT,
    nearest_tm           REAL,
    nearest_identity_fw  REAL,
    nearest_identity_cdr REAL,
    has_experimental_badge INTEGER DEFAULT 0,
    status               TEXT
);

CREATE TABLE IF NOT EXISTS thermostability (
    candidate_id      TEXT PRIMARY KEY REFERENCES candidates(id),
    tm_pred           REAL,
    tm_err            REAL,
    train_nn_distance REAL,
    temstapro         REAL,
    thermompnn_ddg    REAL,
    model_version     TEXT,
    status            TEXT
);

CREATE TABLE IF NOT EXISTS cofolds (
    candidate_id             TEXT NOT NULL REFERENCES candidates(id),
    target                   TEXT NOT NULL REFERENCES targets(target),
    seed                     INTEGER,
    iptm                     REAL,
    ptm                      REAL,
    plddt                    REAL,
    pae                      REAL,
    complex_iptm             REAL,
    prodigy_dg               REAL,
    bsa_total                REAL,
    bsa_per_subunit_json     TEXT,
    n_subunits_contacted     INTEGER,
    epitope_bin              INTEGER,
    interface_residues_json  TEXT,
    path                     TEXT,
    status                   TEXT,
    PRIMARY KEY (candidate_id, target, seed)
);
CREATE INDEX IF NOT EXISTS idx_cofolds_target ON cofolds(target);

CREATE TABLE IF NOT EXISTS scores (
    candidate_id     TEXT PRIMARY KEY REFERENCES candidates(id),
    target           TEXT NOT NULL REFERENCES targets(target),
    axis_affinity    REAL,
    axis_specificity REAL,
    axis_thermo      REAL,
    axis_expression  REAL,
    axis_nativeness  REAL,
    pareto_rank      INTEGER,
    on_front         INTEGER DEFAULT 0,
    desirability     REAL,
    tier_reached     INTEGER,
    axis_affinity_provisional INTEGER DEFAULT 1
);
CREATE INDEX IF NOT EXISTS idx_scores_target ON scores(target);

-- `rank` is the display ordering position within a target. Kept apart from
-- `scores` so that table stays exactly the shape spec Section 7 sets out.
CREATE TABLE IF NOT EXISTS ranks (
    candidate_id TEXT PRIMARY KEY REFERENCES candidates(id),
    target       TEXT,
    rank         INTEGER
);

CREATE TABLE IF NOT EXISTS controls (
    id                     TEXT PRIMARY KEY,
    target                 TEXT,
    name                   TEXT,
    source_doi             TEXT,
    pdb_id                 TEXT,
    sequence               TEXT,
    observed_kd            REAL,
    observed_tm            REAL,
    observed_expression    REAL,
    predicted_affinity     REAL,
    predicted_thermo       REAL,
    predicted_specificity  REAL,
    predicted_expression   REAL,
    kind                   TEXT,   -- known | scrambled | multi_target
    blind                  INTEGER DEFAULT 1,
    assertion              TEXT,
    assertion_passed       INTEGER,
    notes                  TEXT
);

CREATE TABLE IF NOT EXISTS "references" (
    key           TEXT PRIMARY KEY,
    kind          TEXT,
    title         TEXT,
    authors       TEXT,
    year          INTEGER,
    container     TEXT,
    doi           TEXT,
    doi_verified  INTEGER DEFAULT 0,
    url           TEXT,
    repo          TEXT,
    commit_sha    TEXT,
    licence       TEXT,
    version_used  TEXT,
    snapshot_date TEXT,
    role          TEXT,
    sort_order    INTEGER
);

CREATE TABLE IF NOT EXISTS provenance (
    key   TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS qc (
    key    TEXT PRIMARY KEY,
    value  TEXT,
    passed INTEGER
);
"""


def create(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


def main() -> None:
    conn = connect()
    create(conn)
    tables = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    conn.close()
    log(f"schema: {DB_PATH.name} has {len(tables)} tables ({', '.join(tables)})")


if __name__ == "__main__":
    main()
