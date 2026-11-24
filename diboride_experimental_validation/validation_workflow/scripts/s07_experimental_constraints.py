"""Step 7: typed experimental constraint table for ZrB2-NbB2 and ZrB2-TaB2 (item 5).

Every record is assigned exactly one constraint_type, and the table states which part of a thermodynamic model the
record may constrain. The point is to stop melting/solidus data being used as solid-state miscibility boundaries,
and vice versa.

  solidus_incipient_melting   Eckert Pirani incipient melting: lower bound on the solidus at x (liquid+solid model)
  collapse_temperature        Eckert sample collapse: between solidus and liquidus (NOT a liquidus)
  drawn_max_solidus           Eckert's hand-drawn "maximum solidus" curve: author interpretation, not data
  endmember_melting_reference literature end-member melting points (Rudy & Windisch) shown in Eckert's figures
  lattice_parameter           single-set lattice parameters: no boundary information at XRD resolution
  xrd_composition_spread      Otani (10-11) excess broadening -> composition spread (solid-state segregation;
                              T unknown; s11 of the discrepancy analysis)
  near_solidus_solubility     Otani FZ crystal limits (growth from melt): solid solubility near solidus OR growth
                              breakdown; usable only as a hypothesis-dependent bound
  short_anneal_solubility     Vorotilo 1900 C / 10 min: non-equilibrium per authors
  segregation_after_melting   Ness lamellae: segregation exists; equilibrium nature, phases, compositions unknown

Digitization uncertainties are carried unchanged from raw_data/digitized_figures. Eckert's method offset on pure
NbB2 (-57 C collapse, -109 C incipient vs Rudy & Windisch) is recorded as a systematic uncertainty and NOT
subtracted from the data.

Output: results/s07_experimental_constraints.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

WF = Path(__file__).resolve().parents[1]
VAL = WF.parent
DIG = VAL / "raw_data" / "digitized_figures"
CONSTRAINS = {
    "solidus_incipient_melting": ("liquid + AlB2 solid (solidus)", "no"),
    "collapse_temperature": ("liquid + AlB2 solid (between solidus and liquidus)", "no"),
    "drawn_max_solidus": ("none (interpretation)", "no"),
    "endmember_melting_reference": ("end-member liquid/solid lattice stabilities", "no"),
    "lattice_parameter": ("molar volume / Vegard only", "no"),
    "xrd_composition_spread": ("solid-state free energy (gap or segregation), temperature unknown", "indicator only"),
    "near_solidus_solubility": ("solid solvus near solidus IF growth breakdown is thermodynamic", "conditional"),
    "short_anneal_solubility": ("lower bound on solid solubility at 2173 K (non-equilibrium)", "lower bound only"),
    "segregation_after_melting": ("solidification path; mechanism undetermined", "no"),
}


def eckert() -> list[dict]:
    rows = []
    figs = {41: ("ZrNb", "lattice"), 42: ("ZrNb", "melt"), 43: ("ZrTa", "lattice"), 44: ("ZrTa", "melt")}
    for fig, (sid, kind) in figs.items():
        d = pd.read_csv(DIG / f"eckert1966_fig{fig}.csv")
        for r in d.itertuples():
            s = r.series
            if kind == "melt":
                ctype = ("endmember_melting_reference" if s.startswith("Rudy") else "solidus_incipient_melting" if s.startswith("incipient")
                         else "collapse_temperature" if s.startswith("collapse") else "drawn_max_solidus")
                val, unit, unc = r.y_value + 273.15, "K", r.y_uncertainty
            else:
                if "curve" in s:
                    continue  # fitted curve: kept in raw data, not a measurement
                ctype = "lattice_parameter"
                val, unit, unc = r.y_value, f"Angstrom ({s.split()[0]})", r.y_uncertainty
            rows.append({"source": f"Eckert 1966 Fig. {fig}", "point_id": r.point_id, "system": sid, "constraint_type": ctype,
                         "x_solute": r.x_solute, "x_uncertainty": r.x_uncertainty / 100, "value": val, "value_unit": unit,
                         "value_uncertainty": unc,
                         "systematic_note": "Eckert method reads 57-109 C low on pure NbB2 (Fig. 46); not corrected" if kind == "melt" else "",
                         "thermal_history": "hot-pressed, Pirani melted in He; some samples partly unreacted"})
    return rows


def otani() -> list[dict]:
    s11 = pd.read_csv(VAL / "simulation_comparison" / "discrepancy" / "s11_otani_xrd_composition_spread.csv")
    rows = []
    for r in s11[s11.system.isin(["ZrB2-NbB2", "ZrB2-TaB2"])].itertuples():
        rows.append({"source": "Otani 2009 Fig. 1 (via s11)", "point_id": f"OT-{r.system}-{r.x_solute}-{r.cooling[:4]}",
                     "system": "ZrNb" if "Nb" in r.system else "ZrTa", "constraint_type": "xrd_composition_spread",
                     "x_solute": r.x_solute, "x_uncertainty": 0.004, "value": r.dx_equiv_linear, "value_unit": "equivalent dx (linear)",
                     "value_uncertainty": abs(r.dx_equiv_quadrature - r.dx_equiv_linear),
                     "systematic_note": "assumes Cu Ka1 and end-member size/strain width", "thermal_history": f"arc melt, {r.cooling}"})
    return rows


def processed() -> list[dict]:
    rows = []
    sol = pd.read_csv(VAL / "processed_data" / "solubility_limits.csv")
    for r in sol[(sol.host_diboride == "ZrB2") & sol.solute_diboride.isin(["NbB2", "TaB2"])].itertuples():
        ctype = "near_solidus_solubility" if "floating-zone" in str(r.method) else "short_anneal_solubility"
        rows.append({"source": r.source_id, "point_id": r.record_id, "system": "ZrNb" if r.solute_diboride == "NbB2" else "ZrTa",
                     "constraint_type": ctype, "x_solute": r.x_solute_upper_bound, "x_uncertainty": abs(float(r.x_solute_upper_bound) - float(r.x_solute_lower_bound)),
                     "value": r.temperature_K, "value_unit": "K", "value_uncertainty": None, "systematic_note": r.bound_type,
                     "thermal_history": r.temperature_note})
    sep = pd.read_csv(VAL / "processed_data" / "phase_separation.csv")
    for r in sep[sep.system_id.isin(["ZrNb", "ZrTa"])].itertuples():
        if "laser" not in str(r.thermal_history):
            continue
        rows.append({"source": r.source_id, "point_id": r.record_id, "system": r.system_id, "constraint_type": "segregation_after_melting",
                     "x_solute": r.x_groupV, "x_uncertainty": None, "value": r.temperature_K, "value_unit": "K (not measured)",
                     "value_uncertainty": None, "systematic_note": r.equilibrium_or_kinetic_note, "thermal_history": r.thermal_history})
    return rows


def main() -> None:
    df = pd.DataFrame(eckert() + otani() + processed())
    df["constrains_model_part"] = df.constraint_type.map(lambda c: CONSTRAINS[c][0])
    df["usable_as_solid_miscibility_boundary"] = df.constraint_type.map(lambda c: CONSTRAINS[c][1])
    df.to_csv(WF / "results" / "s07_experimental_constraints.csv", index=False)
    print(df.groupby(["system", "constraint_type", "usable_as_solid_miscibility_boundary"]).size().to_string())


if __name__ == "__main__":
    main()
