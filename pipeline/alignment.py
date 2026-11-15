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
# How far the deposited numbering may be shifted from the sequence
# before the placement is refused rather than guessed at.
SHIFT_WINDOW = 40


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


def experimental_structure(path, keys: list[str], placed: str,
                           full: str) -> tuple[str, str]:
    """SSE for a control, from the deposited chain rather than a model.

    A crystal structure is not a model. Residues the density did not resolve
    are simply absent, and the chain carries the entry's own author numbering,
    which is not IMGT and does not reliably start at 1. The deposited chain can
    also hold residues the sequence does not, such as an expression tag: three
    of the ten chains here are longer than the SAbDab-nano sequence they belong
    to.

    So the chain is placed by the one integer shift that makes every residue it
    has agree with the sequence. The shift is searched rather than assumed, the
    search refuses a tie, and a position the chain does not cover is left
    unassigned, which is the honest answer for a disordered loop: that part of
    the structure was not seen.
    """
    if not path or not path.exists():
        return "", "no deposited chain on disk"
    try:
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb

        text = gzip.open(path, "rt").read()
        atoms = pdb.PDBFile.read(io.StringIO(text)).get_structure(model=1)
        sse = struc.annotate_sse(atoms)
        ca = atoms[atoms.atom_name == "CA"]
        if len(sse) != len(ca):
            return "", f"sse {len(sse)} against {len(ca)} residues"
        observed = [(int(residue),
                     struc.info.one_letter_code(name) or "X",
                     SSE_LETTER.get(int(code), "-"))
                    for residue, name, code in zip(ca.res_id, ca.res_name, sse)]
    except Exception as exc:
        return "", f"{type(exc).__name__}"
    if not observed or not full:
        return "", "nothing to place"

    fits: list[tuple[int, int]] = []
    for shift in range(-SHIFT_WINDOW, SHIFT_WINDOW + 1):
        matched = 0
        for residue, aa, _code in observed:
            index = residue - 1 - shift
            if not 0 <= index < len(full):
                continue
            if full[index] != aa:
                matched = -1
                break
            matched += 1
        if matched > 0:
            fits.append((matched, shift))
    if not fits:
        return "", "no placement of the deposited chain agrees with the sequence"
    fits.sort(reverse=True)
    if len(fits) > 1 and fits[0][0] == fits[1][0]:
        return "", "the deposited chain places equally well in two positions"
    covered, shift = fits[0]
    if covered < 0.6 * len(full):
        return "", (f"the deposited chain covers only {covered} of "
                    f"{len(full)} residues")

    offset = len(full) - len(placed)      # residues the numbering could not place
    letters = [UNASSIGNED] * len(placed)
    for residue, _aa, code in observed:
        index = residue - 1 - shift - offset
        if 0 <= index < len(letters):
            letters[index] = code
    if all(letter == UNASSIGNED for letter in letters):
        return "", "the deposited chain covers none of the numbered residues"
    return "".join(letters), (f"the deposited chain, by P-SEA over its CA trace; "
                              f"{covered} of {len(full)} residues resolved")


def reference_chain_path(pdb_id: str, chain: str):
    if not pdb_id or not chain:
        return None
    return STRUCTURES / "_refs" / f"{pdb_id.lower()}_{chain}.pdb.gz"


def model_path(target: str, candidate_id: str):
    path = STRUCTURES / target / f"{candidate_id}.pdb.gz"
    return path if path.exists() else None


def control_rows(conn, target: str) -> list[dict]:
    """The known nanobodies to align beside this target's candidates.

    Two different things, and the panel says which it is showing:

      * a published nanobody with an experimental complex against this target,
        where one exists. EGFR and TNFa have them; GLP1R, PRLR and TL1A have
        none in SAbDab-nano, which is the same coverage gap the Control set
        tab states.
      * failing that, the deposited nanobody chain that this target's own
        shortlist is most often nearest to. It is a framework reference and
        nothing more: it does not bind this target, and the row says so. It is
        chosen by a rule rather than by taste, so it moves when the shortlist
        moves.
    """
    binders = conn.execute(
        "SELECT id, name, pdb_id, chain, sequence, kind FROM controls "
        "WHERE kind = 'known' AND target = ? ORDER BY id", (target,)).fetchall()
    if binders:
        return [dict(r) | {"role": "binder"} for r in binders]

    candidates = conn.execute(
        "SELECT s.nearest_pdb AS pdb_id, s.nearest_chain AS chain, "
        "       COUNT(*) AS n, AVG(s.nearest_identity_fw) AS fw "
        "FROM ranks r JOIN structures s ON s.candidate_id = r.candidate_id "
        "WHERE r.target = ? AND r.rank <= ? AND s.nearest_pdb IS NOT NULL "
        "GROUP BY s.nearest_pdb, s.nearest_chain "
        "ORDER BY n DESC, fw DESC LIMIT 6", (target, shortlist_size())).fetchall()
    best = sequence = sse = None
    for row in candidates:
        path = reference_chain_path(row["pdb_id"], row["chain"])
        found, found_sse = deposited_chain(path)
        if not found:
            continue
        # SAbDab-nano carries conventional VH and TCR variable domains as well
        # as VHHs (see the note in controls.build_known_binders). No candidate
        # in this build has one as its nearest neighbour, but a reference row
        # that is not a nanobody would be worse than none, so it is checked.
        numb = nb.number_many([("q", found)])["q"]
        _intact, n_camelid, _detail = nb.hallmark(numb)
        if n_camelid < 3:
            continue
        best, sequence, sse = row, found, found_sse
        break
    if not sequence:
        log(f"align: {target} has no target-matched control and none of its "
            f"shortlist's commonest experimental neighbours could be read as a "
            f"VHH; the panel gets no control row rather than an invented one")
        return []
    return [{
        "id": f"REF-{target}",
        "name": f"{best['pdb_id']} chain {best['chain']}",
        "pdb_id": best["pdb_id"], "chain": best["chain"],
        "sequence": sequence, "kind": "known", "role": "framework",
        # Already in the chain's own order, so it needs no placing.
        "sse_in_order": sse,
        "n_nearest": best["n"], "fw_identity": round(best["fw"] or 0, 4),
    }]


def deposited_chain(path) -> tuple[str, str]:
    """(observed sequence, one SSE letter per observed residue) for a chain.

    Both are in the order the chain was read, which is what makes them usable
    together: a chain with an unobserved loop yields a sequence that is only
    the residues that were seen, so there is no full sequence to place them
    against and no shift to search for. Position i of one is position i of the
    other, by construction.
    """
    if not path or not path.exists():
        return "", ""
    try:
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb

        atoms = pdb.PDBFile.read(
            io.StringIO(gzip.open(path, "rt").read())).get_structure(model=1)
        sse = struc.annotate_sse(atoms)
        ca = atoms[atoms.atom_name == "CA"]
        if len(sse) != len(ca):
            return "", ""
        return ("".join(struc.info.one_letter_code(n) or "X" for n in ca.res_name),
                "".join(SSE_LETTER.get(int(code), "-") for code in sse))
    except Exception:
        return "", ""


def shortlist_size() -> int:
    from pipeline.common import targets_config

    return int(targets_config().get("run", {}).get("shortlist_per_target", 24))


def main() -> int:
    conn = connect()
    totals = {"rows": 0, "with_sse": 0, "no_model": 0, "mismatch": 0,
              "unnumbered": 0, "unassigned": 0, "controls": 0,
              "controls_with_sse": 0}
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
            elif placed not in sequence:
                # The placed residues are always a contiguous run of the
                # sequence: ANARCI declines residues at the ends, never in the
                # middle. If that ever stops being true the offset below is
                # meaningless, so it is checked rather than assumed.
                raise SystemExit(
                    f"alignment: the residues numbered for {cid} are not a "
                    f"contiguous stretch of its sequence. The grid cannot be "
                    f"written against a sequence it does not match.")
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
            payload.append((cid, "candidate", target, numb.method, len(keys),
                            " ".join(keys), placed, sequence.find(placed),
                            sse, note if sse else "", cid, status))
            totals["rows"] += 1

        # The known controls, through the same numbering as the candidates.
        controls = control_rows(conn, target)
        numbered_controls = nb.number_many(
            [(c["id"], c["sequence"]) for c in controls])
        for control in controls:
            numb = numbered_controls[control["id"]]
            keys = [key_of(pos, ins) for pos, ins, _aa in numb.numbered]
            placed = "".join(aa for _pos, _ins, aa in numb.numbered)
            offset = max(0, control["sequence"].find(placed))
            if control.get("sse_in_order"):
                # A framework reference's sequence is the chain's own observed
                # residues, so the structure needs no placing: the only shift
                # is whatever the numbering declined to place at the front.
                sse = control["sse_in_order"][offset:offset + len(placed)]
                note = ("the deposited chain, by P-SEA over its CA trace; "
                        "every residue of this row is one the structure resolved")
                if len(sse) != len(placed):
                    sse, note = "", "the chain and the numbering disagree in length"
            else:
                sse, note = experimental_structure(
                    reference_chain_path(control.get("pdb_id"),
                                         control.get("chain")),
                    keys, placed, control["sequence"])
            if sse:
                totals["controls_with_sse"] += 1
                totals["unassigned"] += sse.count(UNASSIGNED)
            payload.append((control["id"], "control", target, numb.method,
                            len(keys), " ".join(keys), placed,
                            control["sequence"].find(placed),
                            sse, note if sse else "", control["name"],
                            control["role"]))
            totals["controls"] += 1

        conn.execute("DELETE FROM alignment WHERE target = ?", (target,))
        conn.executemany(
            "INSERT INTO alignment (row_id, kind, target, method, n_residues, "
            "keys, residues, offset, sse, sse_source, label, status) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", payload)
        conn.commit()
        distinct = len({p[5] for p in payload})
        log(f"align: {target} {len(payload)} rows ({len(controls)} of them "
            f"known controls: {', '.join(c['name'] for c in controls) or 'none'}), "
            f"{distinct} distinct IMGT skeletons, "
            f"{sum(1 for p in payload if p[6])} with secondary structure")

    log(f"align: {totals['rows']} rows; secondary structure on "
        f"{totals['with_sse']}; {totals['no_model']} without a model; "
        f"{totals['unassigned']} residues left unassigned; "
        f"{totals['controls']} known control rows, {totals['controls_with_sse']} "
        f"of them with secondary structure from a deposited chain"
        + ("; refused: " + ", ".join(f"{k} x{v}" for k, v in sorted(reasons.items())
                                     if k != "no model") if totals["mismatch"] else ""))
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
