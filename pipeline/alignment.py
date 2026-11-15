"""IMGT-anchored alignment of a target's candidates. Feeds the alignment panel.

Nothing here is aligned. Every candidate in a library is the same fold, and
ANARCI has already placed each of its residues on the IMGT numbering, so two
residues belong in the same column exactly when they carry the same IMGT key.
The alignment is the numbering pivoted, and it is read off rather than
computed.

That is a deliberate choice over a progressive aligner. MUSCLE or MAFFT over
twenty-four VHHs is free to slide a CDR against a framework to buy a cheaper
gap, and the moment it does, a column stops meaning one structural position
across the rows. Reading a conserved framework position off such an alignment
would be reading an artefact of the gap penalty. The IMGT grid cannot do that:
position 42 is position 42 in every row or it is absent.

`tier0` already computes this numbering and keeps only the CDR strings. This
stage keeps the per-residue part, and adds one thing tier0 has no structure
for:

  * per-residue secondary structure from the candidate's own NanoBodyBuilder2
    model, by P-SEA (Labesse 1997) over the CA trace. It is this candidate's
    model rather than a consensus cartoon, which matters most in CDR3, where
    the rows differ and a consensus would be meaningless.

Identities, similarity and the column consensus are not computed here. They
belong to whichever rows the reader has chosen to show, and the reader can
change that, so they are computed in the browser over the current selection.

    pixi run align
"""
from __future__ import annotations

import gzip
import io
import sys

from pipeline import numbering as nb
from pipeline.common import STRUCTURES, TARGETS, connect, log

# Per-residue secondary structure, as one letter. P-SEA resolves helix and
# strand only, so a coil here means "neither", not "measured as random coil".
SSE_LETTER = {-1: "-", 0: "C", 1: "H", 2: "E"}
# A residue the model cannot speak for: no such IMGT key, or the model
# carries a different amino acid there (the germline N-terminal heptad).
UNASSIGNED = "."


def key_of(position: int, insertion: str) -> str:
    """The IMGT key for one residue: '42', or '111A' inside a long CDR3."""
    insertion = (insertion or " ").strip()
    return f"{position}{insertion}"


def secondary_structure(path, keys: list[str], placed: str) -> tuple[str, str]:
    """(one SSE letter per placed residue, a note), or ('', why not).

    The letters are matched to the candidate by IMGT key rather than by
    position, because the model is numbered by the same scheme and the two can
    legitimately disagree about how many residues there are. Fifty-five models
    in this build carry a germline N-terminal heptad in place of a candidate's
    own unnumberable one, which is recorded in `structures.model_note`; those
    residues are a different sequence from the one in the row, so their
    structure is not this candidate's and the column is left unassigned rather
    than filled with a neighbour's answer. A residue takes an SSE letter only
    when the model has that IMGT key and the residue there is the same amino
    acid.
    """
    if not path or not path.exists():
        return "", "no model"
    try:
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb

        text = gzip.open(path, "rt").read() if str(path).endswith(".gz") \
            else path.read_text()
        atoms = pdb.PDBFile.read(io.StringIO(text)).get_structure(model=1)
        sse = struc.annotate_sse(atoms)
        ca = atoms[atoms.atom_name == "CA"]
        if len(sse) != len(ca):
            return "", f"sse {len(sse)} against {len(ca)} residues"
        at_key = {}
        for residue, insertion, name, code in zip(
                ca.res_id, ca.ins_code, ca.res_name, sse):
            at_key[key_of(int(residue), insertion)] = (
                struc.info.one_letter_code(name) or "X",
                SSE_LETTER.get(int(code), "-"))
    except Exception as exc:  # a broken model is not a reason to lose the stage
        return "", f"{type(exc).__name__}"
    letters = []
    for key, aa in zip(keys, placed):
        found = at_key.get(key)
        letters.append(found[1] if found and found[0] == aa else UNASSIGNED)
    if all(letter == UNASSIGNED for letter in letters):
        return "", "no IMGT key in the model matches this row"
    return "".join(letters), "P-SEA over the CA trace of this candidate's model"


def model_path(target: str, candidate_id: str):
    path = STRUCTURES / target / f"{candidate_id}.pdb.gz"
    return path if path.exists() else None


def main() -> int:
    conn = connect()
    totals = {"rows": 0, "with_sse": 0, "no_model": 0, "mismatch": 0,
              "unnumbered": 0, "unassigned": 0}
    reasons: dict[str, int] = {}

    for target in TARGETS:
        rows = conn.execute(
            "SELECT id, sequence FROM candidates WHERE target = ? ORDER BY id",
            (target,)).fetchall()
        numbered = nb.number_many([(r["id"], r["sequence"]) for r in rows])

        payload = []
        for row in rows:
            cid, sequence = row["id"], row["sequence"]
            numb = numbered[cid]
            keys = [key_of(pos, ins) for pos, ins, _aa in numb.numbered]
            placed = "".join(aa for _pos, _ins, aa in numb.numbered)
            status = "ok"
            if not keys:
                status = "unnumbered"
                totals["unnumbered"] += 1
            elif placed != sequence:
                # The numbering may legitimately drop an unnumberable N-terminal
                # extension. Keep the keys, say so, and let the panel show the
                # placed residues only rather than sliding the row by the
                # difference.
                status = f"placed {len(placed)} of {len(sequence)}"

            sse, note = secondary_structure(
                model_path(target, cid), keys, placed)
            totals["unassigned"] += sse.count(UNASSIGNED)
            if sse:
                totals["with_sse"] += 1
            else:
                reasons[note] = reasons.get(note, 0) + 1
                if note == "no model":
                    totals["no_model"] += 1
                else:
                    totals["mismatch"] += 1
            payload.append((cid, target, numb.method, len(keys), " ".join(keys),
                            sse, note if sse else "", status))
            totals["rows"] += 1

        conn.execute("DELETE FROM alignment WHERE target = ?", (target,))
        conn.executemany(
            "INSERT INTO alignment (candidate_id, target, method, n_residues, "
            "keys, sse, sse_source, status) VALUES (?,?,?,?,?,?,?,?)", payload)
        conn.commit()
        distinct = len({p[4] for p in payload})
        log(f"align: {target} {len(payload)} rows, {distinct} distinct IMGT "
            f"skeletons, {sum(1 for p in payload if p[5])} with secondary structure")

    log(f"align: {totals['rows']} rows; secondary structure on "
        f"{totals['with_sse']}; {totals['no_model']} without a model; "
        f"{totals['unassigned']} residues left unassigned because the model "
        f"carries a different amino acid at that IMGT key"
        + ("; refused: " + ", ".join(f"{k} x{v}" for k, v in sorted(reasons.items())
                                     if k != "no model") if totals["mismatch"] else ""))
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
