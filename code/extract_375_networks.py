#!/usr/bin/env python3
"""Export the 375 confirmed-multistable (3,6,3) zero-one networks as (Y, N) matrix pairs.

Reads the pipeline's intermediate index/data files (candidate_indices_*.csv,
potential_multistable_indices.txt, the master enumeration table, odesystem.csv,
check_multistability_*.txt) and writes a single JSON + a human-readable txt file
listing all 375 networks as (Y, N) matrix pairs, with derived reaction equations
and a self-verification pass against the stored ODEs.
"""

import argparse
import json
import re
import sys
from csv import reader as csv_reader
from dataclasses import dataclass, field, asdict
from pathlib import Path

NUM_SPECIES = 3
NUM_REACTIONS = 6
EXPECTED_CANDIDATE_COUNT = 375
SPECIES_NAMES = ("X1", "X2", "X3")
ANOMALOUS_M = 239089  # kept only for reference; all logic is now uniform per-sample-point


class ExtractionError(RuntimeError):
    """Integrity violation that makes it unsafe to proceed."""


@dataclass
class NetworkRecord:
    rank: int
    m: int
    master_index: int
    Y: list
    N: list
    reactions: list = field(default_factory=list)
    ode_check_passed: bool = None
    ode_problems: list = field(default_factory=list)
    maple_verification: dict = None
    multistability_confirmed: bool = None


def load_candidate_m_values(input_dir: Path) -> list:
    paths = sorted(input_dir.glob("candidate_indices_*.csv"))
    if not paths:
        raise ExtractionError(f"no candidate_indices_*.csv files found in {input_dir}")
    values = []
    for p in paths:
        with p.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    values.append(int(line))
    unique_sorted = sorted(set(values))
    if len(values) != len(unique_sorted):
        dupes = sorted({v for v in values if values.count(v) > 1})
        raise ExtractionError(f"duplicate candidate m-values across {paths}: {dupes[:10]}")
    if len(unique_sorted) != EXPECTED_CANDIDATE_COUNT:
        raise ExtractionError(
            f"expected {EXPECTED_CANDIDATE_COUNT} unique candidate m-values, "
            f"found {len(unique_sorted)} across {[p.name for p in paths]}"
        )
    return unique_sorted


def scan_lines_by_number(path: Path, wanted: set, parse) -> dict:
    found = {}
    remaining = set(wanted)
    with path.open(encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            if lineno in remaining:
                found[lineno] = parse(line)
                remaining.discard(lineno)
    if remaining:
        raise ExtractionError(
            f"{path.name}: {len(remaining)} requested line(s) not found, "
            f"e.g. {sorted(remaining)[:10]}"
        )
    return found


def parse_master_row(line: str) -> list:
    ints = [int(tok) for tok in line.strip().split(",")]
    if len(ints) != NUM_SPECIES * NUM_REACTIONS * 2:
        raise ExtractionError(f"expected 36 ints, got {len(ints)}: {line!r}")
    return ints


def parse_ode_row(line: str) -> list:
    row = next(csv_reader([line.strip()]))
    if len(row) != NUM_SPECIES:
        raise ExtractionError(f"expected 3 ODE columns, got {len(row)}: {line!r}")
    return row


def reshape_Y_N(ints: list):
    n = NUM_REACTIONS
    Y = [ints[0:n], ints[n:2 * n], ints[2 * n:3 * n]]
    N = [ints[3 * n:4 * n], ints[4 * n:5 * n], ints[5 * n:6 * n]]
    return Y, N


def format_complex(exponents) -> str:
    bad = [(s, e) for s, e in zip(SPECIES_NAMES, exponents) if e not in (0, 1)]
    if bad:
        raise ExtractionError(
            f"complex exponent outside {{0,1}}: {bad} (indicates a Y/N reshape bug)"
        )
    terms = [s for s, e in zip(SPECIES_NAMES, exponents) if e == 1]
    return " + ".join(terms) if terms else "0"


def derive_reactions(Y: list, N: list) -> list:
    reactions = []
    for j in range(NUM_REACTIONS):
        y_col = [Y[s][j] for s in range(NUM_SPECIES)]
        n_col = [N[s][j] for s in range(NUM_SPECIES)]
        p_col = [y_col[s] + n_col[s] for s in range(NUM_SPECIES)]
        reactant = format_complex(y_col)
        product = format_complex(p_col)
        reactions.append({
            "index": j + 1,
            "rate_constant": f"k{j + 1}",
            "reactant": reactant,
            "product": product,
            "text": f"{reactant} --k{j + 1}--> {product}",
        })
    return reactions


def build_sympy_verifier():
    import sympy as sp

    x_syms = sp.symbols("x1 x2 x3", positive=True)
    k_syms = sp.symbols("k1 k2 k3 k4 k5 k6", positive=True)
    locals_map = {f"x{i + 1}": x_syms[i] for i in range(3)}
    locals_map.update({f"k{i + 1}": k_syms[i] for i in range(6)})

    def reconstruct_ode(Y, N):
        dx = [sp.Integer(0)] * NUM_SPECIES
        for j in range(NUM_REACTIONS):
            monomial = k_syms[j]
            for s in range(NUM_SPECIES):
                if Y[s][j]:
                    monomial *= x_syms[s] ** Y[s][j]
            for s in range(NUM_SPECIES):
                if N[s][j]:
                    dx[s] += N[s][j] * monomial
        return dx

    def verify(Y, N, ode_row):
        reconstructed = reconstruct_ode(Y, N)
        problems = []
        for s in range(NUM_SPECIES):
            stored = sp.sympify(ode_row[s], locals=locals_map)
            diff = sp.expand(reconstructed[s] - stored)
            if diff != 0:
                problems.append(
                    f"dx{s + 1}/dt mismatch: reconstructed={reconstructed[s]}, "
                    f"stored={stored}, diff={diff}"
                )
        return (len(problems) == 0), problems

    return verify


BLOCK_HEADER_RE = re.compile(r"网络索引:\s*(\d+)\s*\|\s*样本点数:\s*(\d+)")
SAMPLE_HEADER_RE = re.compile(r"\[样本点\s*(\d+)\]:\s*\[(.*?)\]")
CONC_STATE_RE = re.compile(r"浓度:\s*\[(.*?)\].*?状态:\s*(\S+)", re.S)
PER_SAMPLE_SUMMARY_RE = re.compile(r"小结:\s*该点正驻点(\d+),\s*稳定\s*(\d+)")
SUMMARY_ROW_RE = re.compile(r"^\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*$")


def parse_kv_list(s: str) -> dict:
    out = {}
    for part in s.split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def parse_check_multistability_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    headers = list(BLOCK_HEADER_RE.finditer(text))
    summary_start = text.find("全局汇总统计表")
    result = {}
    for i, hm in enumerate(headers):
        m_val = int(hm.group(1))
        start = hm.end()
        end = headers[i + 1].start() if i + 1 < len(headers) else (
            summary_start if summary_start != -1 else len(text)
        )
        block = text[start:end]
        sample_headers = list(SAMPLE_HEADER_RE.finditer(block))
        example_parameters, example_steady_states = None, None
        if sample_headers:
            first = sample_headers[0]
            example_parameters = parse_kv_list(first.group(2))
            seg_end = sample_headers[1].start() if len(sample_headers) > 1 else len(block)
            seg = block[first.end():seg_end]
            states = []
            for cm in CONC_STATE_RE.finditer(seg):
                conc = parse_kv_list(cm.group(1))
                try:
                    states.append({
                        "x1": float(conc["x1"]), "x2": float(conc["x2"]), "x3": float(conc["x3"]),
                        "stable": cm.group(2) == "稳定",
                    })
                except (KeyError, ValueError):
                    continue
            example_steady_states = states or None
        per_sample = [(int(p), int(s)) for p, s in PER_SAMPLE_SUMMARY_RE.findall(block)]
        max_stable = max((s for _, s in per_sample), default=0)
        result[m_val] = {
            "source_file": path.name,
            "sample_point_count": None,
            "positive_steadystate_count": None,
            "stable_steadystate_count": None,
            "per_sample_summary": [{"pos": p, "stable": s} for p, s in per_sample],
            "max_stable_per_sample": max_stable,
            "example_parameters": example_parameters,
            "example_steady_states": example_steady_states,
        }
    summary_text = text[summary_start:] if summary_start != -1 else ""
    summary_rows_seen = 0
    for line in summary_text.splitlines():
        sm = SUMMARY_ROW_RE.match(line)
        if not sm:
            continue
        summary_rows_seen += 1
        m_val, spc, psc, ssc = (int(g) for g in sm.groups())
        entry = result.setdefault(m_val, {
            "source_file": path.name, "example_parameters": None, "example_steady_states": None,
        })
        entry["sample_point_count"] = spc
        entry["positive_steadystate_count"] = psc
        entry["stable_steadystate_count"] = ssc
    if summary_rows_seen != len(headers):
        print(
            f"WARNING: {path.name}: {len(headers)} detail block(s) but "
            f"{summary_rows_seen} summary-table row(s) -- counts don't match, "
            f"metadata may be incomplete for this file.",
            file=sys.stderr,
        )
    return result


def load_maple_metadata(input_dir: Path) -> dict:
    merged = {}
    for path in sorted(input_dir.glob("check_multistability_*.txt")):
        merged.update(parse_check_multistability_file(path))
    return merged


def format_matrix_rows(M: list) -> list:
    width = max(len(str(v)) for row in M for v in row)
    return [" ".join(str(v).rjust(width) for v in row) for row in M]


def format_network_block(net: NetworkRecord) -> str:
    lines = [f"Network #{net.rank}  (candidate m={net.m}, master index {net.master_index})"]
    lines.append("  Y (reactant matrix, species x reactions):")
    for row in format_matrix_rows(net.Y):
        lines.append(f"    {row}")
    lines.append("  N (stoichiometric matrix, species x reactions):")
    for row in format_matrix_rows(net.N):
        lines.append(f"    {row}")
    lines.append("  Reactions:")
    reactant_width = max(len(r["reactant"]) for r in net.reactions)
    for r in net.reactions:
        lines.append(f"    {r['reactant'].ljust(reactant_width)} --{r['rate_constant']}--> {r['product']}")
    if not net.ode_check_passed:
        lines.append("  *** WARNING: reconstructed ODE does NOT match odesystem.csv for this network ***")
        for problem in net.ode_problems:
            lines.append(f"      {problem}")
    mv = net.maple_verification
    if mv and mv.get("example_parameters") and mv.get("example_steady_states"):
        spc = mv.get("sample_point_count")
        max_st = mv.get("max_stable_per_sample")
        per_sample = mv.get("per_sample_summary", [])
        sample_detail = " | ".join(f"sample {i+1}: {s['pos']} pos, {s['stable']} stable" for i, s in enumerate(per_sample))
        lines.append(
            f"  Example parameters (Maple sample 1 of {spc}; "
            f"max stable in any single sample: {max_st})"
        )
        lines.append(f"    Per-sample: {sample_detail}")
        params = ", ".join(f"{k}={v}" for k, v in mv["example_parameters"].items())
        lines.append(f"    {params}")
        for i, st in enumerate(mv["example_steady_states"], start=1):
            tag = "stable" if st["stable"] else "unstable"
            lines.append(
                f"    state {i}: x1={st['x1']:.6g}, x2={st['x2']:.6g}, x3={st['x3']:.6g}  [{tag}]"
            )
        if mv.get("meets_stability_threshold_per_sample") is False or net.multistability_confirmed is False:
            lines.append(f"  *** NOTE: multistability_confirmed=False — see audit summary at top of file ***")
    return "\n".join(lines)


def write_json_output(networks: list, out_path: Path) -> None:
    data = []
    for net in networks:
        record = {
            "rank": net.rank,
            "candidate_index": net.m,
            "master_index": net.master_index,
            "Y": net.Y,
            "N": net.N,
            "reactions": [r["text"] for r in net.reactions],
            "ode_check_passed": net.ode_check_passed,
            "multistability_confirmed": net.multistability_confirmed,
            "maple_verification": net.maple_verification,
        }
        data.append(record)
    n_confirmed = sum(1 for net in networks if net.multistability_confirmed)
    n_flagged = len(networks) - n_confirmed
    flagged = [net.m for net in networks if not net.multistability_confirmed]
    output = {
        "_audit_summary": {
            "total_networks_exported": len(networks),
            "multistability_confirmed_count": n_confirmed,
            "multistability_flagged_count": n_flagged,
            "flagged_candidate_indices": flagged,
            "note": (
                "multistability_confirmed is computed per-sample-point (not accumulated across "
                "different parameter values). A network is confirmed multistable iff at least "
                "one sample point (parameter value) has ≥2 stable positive steady states. "
                "The original Maple summary table had a bug that accumulated stable counts "
                "across different sample points, which could be misleading (see m=33449)."
            ),
        },
        "networks": data,
    }
    out_path.write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8")


def write_txt_output(networks: list, out_path: Path, n_passed: int) -> None:
    n_confirmed = sum(1 for net in networks if net.multistability_confirmed)
    blocks = [format_network_block(net) for net in networks]
    flagged = [net.m for net in networks if not net.multistability_confirmed]
    header = (
        f"Audit summary: {n_confirmed}/{len(networks)} networks confirmed multistable "
        f"(>=2 stable steady states in a single sample point).\n"
        f"{len(flagged)} flagged: m={flagged}\n"
        f"Multistability requires >=2 stable positive steady states at the SAME parameter "
        f"value — stable counts from different sample points must not be added together."
    )
    footer = f"{n_passed}/{len(networks)} passed ODE self-check"
    out_path.write_text(header + "\n\n" + "\n\n".join(blocks) + "\n\n" + footer + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("/Users/zhangjiandong/Downloads/multistability"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--skip-verification", action="store_true")
    parser.add_argument("--skip-metadata", action="store_true")
    args = parser.parse_args(argv)

    print(f"Reading pipeline data from {args.input_dir}")
    m_values = load_candidate_m_values(args.input_dir)
    print(f"Loaded {len(m_values)} unique candidate m-values.")

    m_to_master = scan_lines_by_number(
        args.input_dir / "potential_multistable_indices.txt", set(m_values), lambda l: int(l.strip())
    )
    master_indices_needed = set(m_to_master.values())
    master_to_ints = scan_lines_by_number(
        args.input_dir / "363_2-algorithm1-4.txt", master_indices_needed, parse_master_row
    )

    networks = []
    for rank, m in enumerate(m_values, start=1):
        master_index = m_to_master[m]
        Y, N = reshape_Y_N(master_to_ints[master_index])
        reactions = derive_reactions(Y, N)
        networks.append(NetworkRecord(rank=rank, m=m, master_index=master_index, Y=Y, N=N, reactions=reactions))

    # duplicate master_index check (warn only, never drop data)
    seen_master = {}
    for net in networks:
        seen_master.setdefault(net.master_index, []).append(net.m)
    dupes = {mi: ms for mi, ms in seen_master.items() if len(ms) > 1}
    if dupes:
        print(f"WARNING: {len(dupes)} master_index value(s) shared by multiple candidates: {dupes}", file=sys.stderr)

    n_passed = len(networks)
    if args.skip_verification:
        print("Skipping ODE self-verification (--skip-verification).")
        for net in networks:
            net.ode_check_passed = None
    else:
        print("Verifying reconstructed ODEs against odesystem.csv with sympy...")
        m_to_ode = scan_lines_by_number(args.input_dir / "odesystem.csv", set(m_values), parse_ode_row)
        verify = build_sympy_verifier()
        n_passed = 0
        for net in networks:
            passed, problems = verify(net.Y, net.N, m_to_ode[net.m])
            net.ode_check_passed = passed
            net.ode_problems = problems
            if passed:
                n_passed += 1
            else:
                print(f"MISMATCH: rank {net.rank} (m={net.m}, master={net.master_index}):", file=sys.stderr)
                for p in problems:
                    print(f"    {p}", file=sys.stderr)
        print(f"ODE self-check: {n_passed}/{len(networks)} passed.")

    if args.skip_metadata:
        print("Skipping Maple metadata enrichment (--skip-metadata).")
        for net in networks:
            net.multistability_confirmed = None
    else:
        try:
            metadata = load_maple_metadata(args.input_dir)
            for net in networks:
                mv = metadata.get(net.m)
                if mv is None:
                    continue
                mv = dict(mv)
                # Compute multistability_confirmed from per-sample-point data
                # (NOT from summary table, which had the accumulation bug)
                max_stable = mv.get("max_stable_per_sample", 0)
                net.multistability_confirmed = (max_stable >= 2)
                # Renamed from meets_stability_threshold to be explicit about per-sample nature
                mv["meets_stability_threshold_per_sample"] = (max_stable >= 2)
                net.maple_verification = mv
        except Exception as exc:  # noqa: BLE001 - bonus metadata must never block core output
            print(f"WARNING: Maple metadata parsing failed ({exc}); proceeding without it.", file=sys.stderr)
            for net in networks:
                net.maple_verification = None
                net.multistability_confirmed = None

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "375_multistable_networks.json"
    txt_path = args.output_dir / "375_multistable_networks.txt"

    n_confirmed = sum(1 for net in networks if net.multistability_confirmed)
    print(f"Multistability audit: {n_confirmed}/{len(networks)} confirmed (per-sample-point check).")
    if n_confirmed != len(networks):
        flagged = [net.m for net in networks if not net.multistability_confirmed]
        print(f"  Flagged (no single sample point has >=2 stable): {flagged}")

    write_json_output(networks, json_path)
    write_txt_output(networks, txt_path, n_passed)
    print(f"Wrote {json_path}")
    print(f"Wrote {txt_path}")

    if not args.skip_verification and n_passed != len(networks):
        print(f"FAILED: only {n_passed}/{len(networks)} networks passed ODE self-check.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
