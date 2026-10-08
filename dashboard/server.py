#!/usr/bin/env python3
"""BLADE Dashboard Server.

Serves a live dashboard for monitoring BLADE pipeline runs.
Reads the same TOML config as full_framework.py — no base code modified.

Usage:
    python dashboard/server.py examples/full_framework.toml
    python dashboard/server.py examples/full_framework.toml --port 8080 --log nohup.out
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import sys
import time
import tomllib
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def _load_config(paths: list[Path]) -> dict:
    cfg: dict = {}
    for p in paths:
        with p.open("rb") as f:
            extra = tomllib.load(f)
        cfg = _deep_merge(cfg, extra)
    return cfg


def _deep_merge(base: dict, override: dict) -> dict:
    result = dict(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = val
    return result


def _resolve_paths(cfg: dict, toml_path: Path) -> tuple[Path, Path]:
    """Return (blade_root, files_dir) from config."""
    blade_root = Path(cfg.get("paths", {}).get("blade_root", toml_path.parent.parent)).expanduser().resolve()
    files_dir = Path(cfg.get("paths", {}).get("files_dir", blade_root / "Files")).expanduser().resolve()
    return blade_root, files_dir


# ---------------------------------------------------------------------------
# ATAT bestsqs.out → XYZ parser
# ---------------------------------------------------------------------------

def _build_elem_map(rndstr_path: Path) -> dict[str, str]:
    """Parse rndstr.in → ATAT component label → real element.

    Returns e.g. {'a_A': 'Hf', 'a_B': 'Cr'}.  Pure sites (no '=') are skipped
    because they appear as the element symbol directly in bestsqs.out.
    """
    try:
        text = rndstr_path.read_text(errors="ignore")
    except OSError:
        return {}
    lines = [l.strip() for l in text.strip().splitlines()
             if l.strip() and not l.startswith("#")]
    if len(lines) < 5:
        return {}
    mapping: dict[str, str] = {}
    sl_letters = "abcdefghijklmnop"
    sl_idx = 0
    for line in lines[4:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            float(parts[0])
        except ValueError:
            continue
        spec = parts[3]
        if "=" in spec:
            sl = sl_letters[sl_idx] if sl_idx < len(sl_letters) else f"s{sl_idx}"
            for ci, token in enumerate(spec.split(",")):
                elem = token.split("=")[0].strip()
                mapping[f"{sl}_{chr(ord('A') + ci)}"] = elem
            sl_idx += 1
    return mapping


# One element per letter A-Z, in periodic-table order for visual variety.
# Used when rndstr.in is unavailable or missing a component.
_ALPHA_ELEMS = [
    "Al", "B",  "Ca", "Cr", "Cu", "Fe", "Ga", "Hf", "In", "Ir",
    "K",  "La", "Mn", "N",  "O",  "P",  "Pt", "Re", "S",  "Ti",
    "U",  "V",  "W",  "Xe", "Y",  "Zr",
]  # A→Al, B→B, C→Ca, D→Cr, E→Cu, F→Fe, G→Ga, H→Hf, ...


def _atat_species(raw: str, elem_map: dict[str, str] | None = None) -> str:
    """Map ATAT alloy notation 'a_A' → real element (rndstr.in map, else alphabet fallback)."""
    if elem_map and raw in elem_map:
        return elem_map[raw]
    m = re.match(r"^[a-z]_([A-Z])$", raw)
    if m:
        idx = ord(m.group(1)) - ord("A")
        return _ALPHA_ELEMS[idx] if idx < len(_ALPHA_ELEMS) else "Fe"
    return raw


def parse_atat_structure(text: str, elem_map: dict[str, str] | None = None) -> str | None:
    """Convert ATAT str.out / bestsqs.out text to XYZ format string."""
    lines = [l.strip() for l in text.strip().splitlines() if l.strip() and not l.startswith("#")]
    if len(lines) < 7:
        return None
    try:
        # bestsqs.out format: lines 0-2 = parent lattice, lines 3-5 = supercell vectors.
        # Atom fractional coords are in the PARENT lattice basis (lines 0-2).
        lat = []
        for i in range(3):
            lat.append([float(x) for x in lines[i].split()[:3]])

        # Parse atoms (remaining lines: x y z species)
        atoms: list[tuple[str, list[float]]] = []
        for line in lines[6:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            # Some ATAT files: x y z species; others: species x y z
            try:
                frac = [float(parts[0]), float(parts[1]), float(parts[2])]
                species = _atat_species(parts[3], elem_map)
            except ValueError:
                # species first format
                species = _atat_species(parts[0], elem_map)
                frac = [float(parts[1]), float(parts[2]), float(parts[3])]

            # Convert fractional → Cartesian
            cart = [
                frac[0] * lat[0][0] + frac[1] * lat[1][0] + frac[2] * lat[2][0],
                frac[0] * lat[0][1] + frac[1] * lat[1][1] + frac[2] * lat[2][1],
                frac[0] * lat[0][2] + frac[1] * lat[1][2] + frac[2] * lat[2][2],
            ]
            atoms.append((species, cart))

        if not atoms:
            return None

        xyz_lines = [str(len(atoms)), "BLADE dashboard"]
        for species, (x, y, z) in atoms:
            xyz_lines.append(f"{species} {x:.6f} {y:.6f} {z:.6f}")
        return "\n".join(xyz_lines)
    except Exception:
        return None


def parse_bestcorr(text: str) -> float | None:
    """Extract objective value from bestcorr*.out."""
    for line in text.strip().splitlines():
        parts = line.strip().split()
        if parts:
            try:
                return float(parts[-1])
            except ValueError:
                continue
    return None


# ---------------------------------------------------------------------------
# Filesystem scanners
# ---------------------------------------------------------------------------

def scan_sqs(files_dir: Path) -> list[dict]:
    """Scan Files/SQS/ for active mcsqs runs.

    mcsqs writes bestcorr<N>.out / bestsqs<N>.out into the directory it is
    spawned from, which may be a per-composition subdir (sqsdb_lev=*) rather
    than the top-level lattice dir.  We collect runs from any depth.
    """
    sqs_root = files_dir / "SQS"
    results = []
    if not sqs_root.exists():
        return results

    for sqs_dir in sorted(sqs_root.iterdir()):
        if not sqs_dir.is_dir():
            continue

        # Group bestcorr files by the directory they live in (could be subdir)
        by_parent: dict[Path, list] = {}
        for bestcorr in sorted(sqs_dir.rglob("bestcorr*.out")):
            by_parent.setdefault(bestcorr.parent, []).append(bestcorr)

        # Flatten: collect all runs across all subdirs
        runs = []
        for parent_dir, corr_files in sorted(by_parent.items()):
            rel = parent_dir.relative_to(sqs_dir)
            for bestcorr in sorted(corr_files):
                ip_str = re.sub(r"\D", "", bestcorr.stem.replace("bestcorr", ""))
                ip = int(ip_str) if ip_str else 0
                obj = None
                try:
                    obj = parse_bestcorr(bestcorr.read_text(errors="ignore"))
                except OSError:
                    pass
                # mcsqs with no -ip writes "bestsqs.out"; with -ip=N writes "bestsqs1.out"..."bestsqsN.out"
                bestsqs = parent_dir / f"bestsqs{ip_str}.out"
                runs.append({
                    "ip": ip,
                    "subdir": str(rel),
                    "bestcorr_path": str(bestcorr),
                    "bestsqs_path": str(bestsqs) if bestsqs.exists() else None,
                    "objective": obj,
                    "mtime": bestcorr.stat().st_mtime,
                })

        final = sqs_dir / "bestsqs.out"
        done = final.exists()

        # Active composition: the subdir with the most-recently-modified bestcorr
        active_comp = None
        if runs:
            latest = max(runs, key=lambda r: r["mtime"])
            active_comp = latest["subdir"]

        results.append({
            "name": sqs_dir.name,
            "path": str(sqs_dir),
            "runs": runs,
            "done": done,
            "final_path": str(final) if done else None,
            "active_comp": active_comp,
        })
    return results


def scan_comps(files_dir: Path, comps_folder: str = "Comps") -> list[dict]:
    """Scan Files/Comps/ for system status."""
    comps_root = files_dir / comps_folder
    results = []
    if not comps_root.exists():
        return results

    for comp_dir in sorted(comps_root.iterdir()):
        if not comp_dir.is_dir():
            continue
        tdb_files = list(comp_dir.glob("*.tdb"))
        energy_files = list(comp_dir.rglob("energy"))
        traj_files = list(comp_dir.rglob("*.xyz"))
        contcar_files = list(comp_dir.rglob("CONTCAR"))

        results.append({
            "name": comp_dir.name,
            "path": str(comp_dir),
            "tdb_done": bool(tdb_files),
            "tdb_files": [str(f) for f in tdb_files],
            "energy_count": len(energy_files),
            "has_trajectory": bool(traj_files),
            "trajectory_paths": [str(f) for f in traj_files],
            "contcar_paths": [str(f) for f in contcar_files],
            "mtime": comp_dir.stat().st_mtime,
        })
    return results


def scan_plots(files_dir: Path) -> list[dict]:
    """Scan for all PNG/GIF outputs."""
    plots = []
    for ext in ("*.png", "*.gif"):
        for f in sorted(files_dir.rglob(ext)):
            plots.append({
                "name": f.name,
                "path": str(f),
                "rel": str(f.relative_to(files_dir)),
                "mtime": f.stat().st_mtime,
                "size": f.stat().st_size,
            })
    plots.sort(key=lambda x: x["mtime"], reverse=True)
    return plots


def scan_compositions(system_dir: Path) -> list[dict]:
    """List individual SQS composition subdirs within a system (each has a CONTCAR)."""
    results = []
    for contcar in sorted(system_dir.rglob("CONTCAR")):
        parent = contcar.parent
        traj = sorted(parent.rglob("*.xyz"))
        energy_file = parent / "energy"
        energy = None
        if energy_file.exists():
            try:
                energy = float(energy_file.read_text().strip().split()[0])
            except (ValueError, IndexError):
                pass
        results.append({
            "name": parent.name,
            "path": str(parent),
            "contcar_path": str(contcar),
            "has_trajectory": bool(traj),
            "trajectory_paths": [str(f) for f in traj],
            "energy": energy,
        })
    return results


def scan_energy(comp_path: Path) -> list[dict]:
    """Parse energy files from a composition directory."""
    results = []
    for energy_file in sorted(comp_path.rglob("energy")):
        try:
            val = float(energy_file.read_text().strip().split()[0])
            results.append({"path": str(energy_file), "energy": val, "rel": str(energy_file.relative_to(comp_path))})
        except (OSError, ValueError, IndexError):
            pass
    return results


def infer_stage_status(files_dir: Path, cfg: dict) -> dict:
    """Infer pipeline stage completion from filesystem."""
    comps_folder = cfg.get("tdb", {}).get("comps_folder", "Comps")
    sqs_root = files_dir / "SQS"
    comps_root = files_dir / comps_folder
    phase_diag_root = files_dir / cfg.get("phase_plots", {}).get("output_folder", "Phase_Diagrams")

    comps = scan_comps(files_dir, comps_folder)
    n_done = sum(1 for c in comps if c["tdb_done"])
    n_total = len(comps)

    return {
        "sqs_generation": {
            "started": sqs_root.exists(),
            "dirs": len(list(sqs_root.iterdir())) if sqs_root.exists() else 0,
        },
        "tdb_fitting": {
            "started": comps_root.exists(),
            "systems_total": n_total,
            "systems_done": n_done,
            "all_done": n_total > 0 and n_done == n_total,
        },
        "phase_diagrams": {
            "started": phase_diag_root.exists(),
            "plot_count": len(list(phase_diag_root.rglob("*.png"))) if phase_diag_root.exists() else 0,
        },
        "composition_list": (files_dir / "composition_list.xlsx").exists(),
    }


def read_log_tail(log_path: Path, n_lines: int = 80) -> list[str]:
    if not log_path or not log_path.exists():
        return []
    try:
        text = log_path.read_text(errors="ignore")
        return text.splitlines()[-n_lines:]
    except OSError:
        return []


# ---------------------------------------------------------------------------
# HTTP handler
# ---------------------------------------------------------------------------

class DashboardHandler(BaseHTTPRequestHandler):
    files_dir: Path
    blade_root: Path
    cfg: dict
    log_path: Path | None

    def log_message(self, fmt, *args):
        pass  # suppress request logging

    def _json(self, data: Any, status: int = 200) -> None:
        body = json.dumps(data, default=str).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(body))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def _serve_file(self, path: Path, content_type: str = "text/plain") -> None:
        try:
            data = path.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", len(data))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(data)
        except OSError:
            self._json({"error": "not found"}, 404)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        qs = parse_qs(parsed.query)
        path = parsed.path.rstrip("/") or "/"

        if path == "/" or path == "/index.html":
            html = Path(__file__).parent / "index.html"
            self._serve_file(html, "text/html; charset=utf-8")

        elif path == "/api/status":
            self._json({
                "stages": infer_stage_status(self.files_dir, self.cfg),
                "sqs": scan_sqs(self.files_dir),
                "systems": scan_comps(self.files_dir, self.cfg.get("tdb", {}).get("comps_folder", "Comps")),
                "timestamp": time.time(),
            })

        elif path == "/api/sqs":
            self._json(scan_sqs(self.files_dir))

        elif path == "/api/structure":
            file_path = qs.get("path", [None])[0]
            if not file_path:
                self._json({"error": "missing path"}, 400)
                return
            p = Path(file_path)
            if not p.exists():
                self._json({"error": "not found"}, 404)
                return
            elem_map = _build_elem_map(p.parent / "rndstr.in")
            xyz = parse_atat_structure(p.read_text(errors="ignore"), elem_map)
            self._json({"xyz": xyz, "path": str(p)})

        elif path == "/api/trajectory":
            file_path = qs.get("path", [None])[0]
            if not file_path:
                self._json({"error": "missing path"}, 400)
                return
            p = Path(file_path)
            if not p.exists():
                self._json({"error": "not found"}, 404)
                return
            self._json({"xyz": p.read_text(errors="ignore"), "path": str(p)})

        elif path == "/api/compositions":
            dir_path = qs.get("dir", [None])[0]
            if not dir_path:
                self._json({"error": "missing dir"}, 400)
                return
            self._json(scan_compositions(Path(dir_path)))

        elif path == "/api/energy":
            dir_path = qs.get("dir", [None])[0]
            if not dir_path:
                self._json({"error": "missing dir"}, 400)
                return
            self._json(scan_energy(Path(dir_path)))

        elif path == "/api/plots":
            self._json(scan_plots(self.files_dir))

        elif path == "/api/image":
            file_path = qs.get("path", [None])[0]
            if not file_path:
                self._json({"error": "missing path"}, 400)
                return
            p = Path(file_path)
            if not p.exists():
                self._json({"error": "not found"}, 404)
                return
            ext = p.suffix.lower()
            mime = "image/gif" if ext == ".gif" else "image/png"
            b64 = base64.b64encode(p.read_bytes()).decode()
            self._json({"data": b64, "mime": mime, "name": p.name})

        elif path == "/api/log":
            lines = read_log_tail(self.log_path)
            self._json({"lines": lines})

        elif path == "/api/config":
            self._json({
                "files_dir": str(self.files_dir),
                "blade_root": str(self.blade_root),
                "elements": self.cfg.get("elements", {}),
                "tdb": {k: v for k, v in self.cfg.get("tdb", {}).items() if k != "primary_elements"},
                "stages": self.cfg.get("stages", {}),
            })

        else:
            self._json({"error": "not found"}, 404)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="BLADE dashboard server")
    parser.add_argument("config", nargs="+", type=Path, help="TOML config file(s), same as full_framework.py")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--log", type=Path, default=None, help="Path to nohup.out log file")
    args = parser.parse_args()

    cfg = _load_config(args.config)
    blade_root, files_dir = _resolve_paths(cfg, args.config[0])

    print(f"BLADE Dashboard")
    print(f"  blade_root : {blade_root}")
    print(f"  files_dir  : {files_dir}")
    print(f"  log        : {args.log or '(none)'}")
    print(f"  http://localhost:{args.port}")
    print()

    DashboardHandler.files_dir = files_dir
    DashboardHandler.blade_root = blade_root
    DashboardHandler.cfg = cfg
    DashboardHandler.log_path = args.log

    server = HTTPServer(("", args.port), DashboardHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")


if __name__ == "__main__":
    main()
