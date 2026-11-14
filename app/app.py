"""WHITTLE web app. Spec Section 8; DESIGN.md is authoritative for layout.

Flask plus Jinja, read-only SQLite, no request-time compute. Every chart and
every table reads a static JSON bundle written by `pipeline/build_db.py` on the
Studio. SQLite is touched for one thing only: the single candidate row the
drawer opens with.

Structure lifted from debris/app/app.py and GOBSMACKED/app (see REUSE.md); the
DEBRIS design system is replaced by DESIGN.md's.
"""
from __future__ import annotations

import csv
import io
import json
import os
import sqlite3
import tomllib
from pathlib import Path

from flask import (Flask, Response, abort, g, jsonify, render_template, request,
                   send_from_directory, url_for)

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = Path(os.environ.get("WHITTLE_DB", ROOT / "data" / "whittle.sqlite"))
STATIC_DATA = Path(os.environ.get(
    "WHITTLE_STATIC_DATA", Path(__file__).resolve().parent / "static" / "data"))
STRUCTURES = Path(os.environ.get("WHITTLE_STRUCTURES", ROOT / "data" / "structures"))

# Public but unlisted while the data is provisional (spec 8.4). Phase 4's final
# deploy sets WHITTLE_INDEXABLE=1 and both the robots rule and the noindex meta
# tag come off together.
INDEXABLE = os.environ.get("WHITTLE_INDEXABLE", "0") == "1"

TABS = [("Candidates", "/"), ("Control set", "/controls"), ("About", "/about")]

BACKRONYM = [
    ("W", "hole-repertoire"), ("H", "ierarchical"), ("I", "n-silico"),
    ("T", "riage:"), ("T", "hermostability,"), ("L", "iabilities,"),
    ("E", "xpression"),
]

app = Flask(__name__)
app.config["JSON_SORT_KEYS"] = False


# ----------------------------------------------------------------------- plumbing

def db() -> sqlite3.Connection:
    if "db" not in g:
        if not DB_PATH.exists():
            abort(503, "database not present")
        conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        g.db = conn
    return g.db


@app.teardown_appcontext
def close_db(_exc) -> None:
    if conn := g.pop("db", None):
        conn.close()


def bundle(name: str) -> dict:
    """A static bundle, read fresh each request so a data swap needs no restart."""
    path = STATIC_DATA / f"{name}.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


def meta() -> dict:
    data = bundle("meta")
    if not data:
        # The app must still render before the first build-db run, so a visitor
        # sees the shell and a plain "data pending" line rather than a 500.
        data = {"run_id": "pending", "built_date": "", "n_candidates": 0,
                "n_targets": 0, "targets": [], "columns": {"total": 0, "groups": [],
                "columns": [], "defaults": []}, "axes": {}, "axes_pending": [],
                "tier_counts": {}}
    return data


def asset_version(relative: str) -> str:
    """`?v=<mtime>` cache busting (spec 8.1, DESIGN 3.6)."""
    path = Path(app.static_folder) / relative
    try:
        return str(int(path.stat().st_mtime))
    except OSError:
        return "0"


@app.context_processor
def inject() -> dict:
    return {
        "tabs": TABS,
        "active": request.path,
        "meta": meta(),
        "backronym": BACKRONYM,
        "asset_version": asset_version,
        "indexable": INDEXABLE,
    }


# ------------------------------------------------------------------------- pages

@app.route("/")
def candidates():
    return render_template("candidates.html", title="Candidates")


@app.route("/controls")
def controls():
    return render_template("controls.html", title="Control set",
                           controls=bundle("controls"), qc=bundle("qc"))


@app.route("/about")
def about():
    data = bundle("about")
    diagram = ROOT / "app" / "static" / "img" / "workflow.svg"
    return render_template(
        "about.html", title="About and methods", about=data, qc=bundle("qc"),
        axes=bundle("axes"),
        diagram_svg=diagram.read_text() if diagram.exists() else "")


# --------------------------------------------------------------------------- api

@app.route("/api/candidate/<cid>.json")
def api_candidate(cid: str):
    """One candidate, everything about it. The drawer's single fetch on open."""
    conn = db()
    row = conn.execute("""
        SELECT c.*, r.rank, s.axis_affinity, s.axis_specificity, s.axis_thermo,
               s.axis_expression, s.axis_nativeness, s.pareto_rank, s.on_front,
               s.desirability, s.tier_reached, s.axis_affinity_provisional,
               st.path AS structure_path, st.refined, st.mean_error, st.cdr3_error,
               st.psh, st.ppc, st.pnc, st.tnp_flags, st.tnp_version, st.swi,
               st.protein_sol, st.aggrescan3d, st.hydrophobic_patch_area,
               st.charged_patch_area, st.sasa_total, st.nativeness_ppl,
               st.oasis_humanness, st.nearest_pdb, st.nearest_tm,
               st.nearest_identity_fw, st.nearest_identity_cdr,
               st.nearest_chain, st.nearest_rmsd, st.superposed,
               st.reference_file,
               st.has_experimental_badge, st.model_note,
               t.tm_pred, t.tm_err, t.train_nn_distance, t.temstapro,
               t.thermompnn_ddg, t.model_version
        FROM candidates c
        LEFT JOIN ranks r ON r.candidate_id = c.id
        LEFT JOIN scores s ON s.candidate_id = c.id
        LEFT JOIN structures st ON st.candidate_id = c.id
        LEFT JOIN thermostability t ON t.candidate_id = c.id
        WHERE c.id = ?""", (cid,)).fetchone()
    if row is None:
        abort(404)
    rec = dict(row)
    rec["liabilities_detail"] = [dict(r) for r in conn.execute(
        "SELECT kind, region, position, detail FROM liabilities "
        "WHERE candidate_id = ? ORDER BY position", (cid,))]
    rec["cofolds"] = [dict(r) for r in conn.execute(
        "SELECT seed, iptm, ptm, plddt, pae, prodigy_dg, bsa_total, "
        "bsa_per_subunit_json, n_subunits_contacted, epitope_bin, "
        "interface_residues_json, status FROM cofolds "
        "WHERE candidate_id = ? ORDER BY seed", (cid,))]
    target = conn.execute("SELECT * FROM targets WHERE target = ?",
                          (rec["target"],)).fetchone()
    rec["target_detail"] = dict(target) if target else {}
    # Structure URL only when the file is actually on this host.
    rec["structure_url"] = None
    rec["superposed_url"] = None
    if rec.get("structure_path"):
        # The .pdb names, not the .pdb.gz ones: see the `structure` route.
        if (STRUCTURES / rec["target"] / f"{cid}.pdb.gz").exists():
            rec["structure_url"] = url_for("structure", target=rec["target"],
                                           name=f"{cid}.pdb")
        # The model already moved onto its nearest experimental neighbour's
        # frame by TM-align at build time, so the drawer loads two files and
        # does no geometry of its own.
        if (STRUCTURES / rec["target"] / f"{cid}.sup.pdb.gz").exists():
            rec["superposed_url"] = url_for("structure", target=rec["target"],
                                            name=f"{cid}.sup.pdb")
    # The co-folded complex: the nanobody docked on the real antigen construct,
    # which is the thing Tier 2 exists to show.
    rec["complex_url"] = None
    if (STRUCTURES / rec["target"] / f"{cid}.complex.pdb.gz").exists():
        rec["complex_url"] = url_for("structure", target=rec["target"],
                                     name=f"{cid}.complex.pdb")
    # The neighbour chain on its own, shared by every candidate that matches
    # it, so the viewer loads two local files and never calls out to RCSB.
    rec["reference_url"] = None
    if rec.get("reference_file"):
        if (STRUCTURES / "_refs" / rec["reference_file"]).exists():
            rec["reference_url"] = url_for(
                "structure", target="_refs",
                name=rec["reference_file"].replace(".pdb.gz", ".pdb"))
    return jsonify(rec)


@app.route("/structures/<target>/<name>")
def structure(target: str, name: str):
    """Predicted models, stored gzipped and requested as plain `.pdb`.

    The URL asks for `<id>.pdb` and the file on disk is `<id>.pdb.gz`. That is
    deliberate: in production nginx's `gzip_static always` picks up the `.gz`
    for a request to the uncompressed name and sets `Content-Encoding: gzip`
    itself, so the browser decompresses transparently and Mol* reads plain PDB.
    Serving `<id>.pdb.gz` at its own URL instead would hand the browser raw
    gzip bytes with no encoding header and the viewer would see rubbish. This
    route is the development equivalent and sets the same headers.
    """
    if "/" in name or ".." in name or not name.endswith(".pdb"):
        # `<id>.sup.pdb` is the superposed copy and ends in .pdb like any other.
        abort(404)
    directory = STRUCTURES / target
    gzipped = directory / f"{name}.gz"
    if not gzipped.exists():
        abort(404)
    response = send_from_directory(directory, f"{name}.gz")
    response.headers["Content-Encoding"] = "gzip"
    response.headers["Content-Type"] = "chemical/x-pdb"
    response.headers["Content-Disposition"] = f'inline; filename="{name}"'
    return response


@app.route("/api/export.csv")
def api_export():
    """CSV of the current view. `target` is required: there is no pooled export."""
    target = request.args.get("target", "")
    valid = {t["target"] for t in meta().get("targets", [])}
    if target not in valid:
        abort(400, "a known target is required: there is no pooled view")
    data = bundle(f"candidates_{target}")
    rows = data.get("rows", [])
    wanted = [k for k in (request.args.get("columns", "").split(",")) if k]
    ids = {i for i in (request.args.get("ids", "").split(",")) if i}
    if ids:
        rows = [r for r in rows if r["id"] in ids]
    keys = wanted or (list(rows[0]) if rows else [])
    out = io.StringIO()
    out.write(f"# WHITTLE {meta().get('run_id')} - {target} - percentiles are "
              f"within this target and are not comparable between targets\n")
    writer = csv.DictWriter(out, fieldnames=keys, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return Response(out.getvalue(), mimetype="text/csv", headers={
        "Content-Disposition": f'attachment; filename="whittle_{target}.csv"'})


@app.route("/healthz")
def healthz():
    data = meta()
    return jsonify({"ok": True, "run_id": data.get("run_id"),
                    "built_at": data.get("built_at"),
                    "git_sha": data.get("git_sha"),
                    "n_candidates": data.get("n_candidates"),
                    "db": DB_PATH.exists()})


@app.route("/robots.txt")
def robots():
    if INDEXABLE:
        body = "User-agent: *\nAllow: /\n"
    else:
        # Public but unlisted while the data is provisional (spec 8.4).
        body = "User-agent: *\nDisallow: /\n"
    return Response(body, mimetype="text/plain")


@app.errorhandler(404)
def not_found(_e):
    return render_template("error.html", title="Not found", code=404,
                           message="No such page or candidate."), 404


@app.errorhandler(503)
def unavailable(e):
    return render_template("error.html", title="Data pending", code=503,
                           message=str(getattr(e, "description", e))), 503


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("PORT", "8021")), debug=True)
