"""MOZART-TS1: the mechanism, its rate constants, the cells, and the scenarios.

The mechanism is read from data/ts1.json, a MICM-style configuration of TS1
(MICM's test/benchmark/ts1_mechanism.hpp at commit e15898ce, rate constants in
mol m-3 units). Rate constants follow the baseline's formulas exactly:

    Arrhenius   k = A exp(C/T) (T/D)^B (1 + E p)
    Troe        k0 = k0_A exp(k0_C/T) (T/300)^k0_B, kinf likewise, r = k0 [M]/kinf,
                k = k0 [M] / (1 + r) * Fc^(N / (N + log10(r)^2))
    user        k = (rate parameter) * (scaling factor)
    surface     k = 4 N pi r^2 / (r / D + 4 / (v gamma)),  v = sqrt(8 R T / (pi Mw))

with [M] = p / (R T). Third bodies (M) are folded into k. A scenario sets every
photolysis parameter (label "j...") to its cell value times f(t) at the start of
each interval (benchmarks/tier3/scenarios/README.md).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from .chem import ChemSystem, Mechanism
from .counts import Counts

DATA = Path(__file__).with_name("data") / "ts1.json"
BOLTZMANN = 1.380649e-23
AVOGADRO = 6.02214076e23
R_GAS = BOLTZMANN * AVOGADRO


class TS1:
    def __init__(self, path: Path = DATA):
        d = json.loads(Path(path).read_text())
        cd = d["camp-data"]
        specs = [e for e in cd if e.get("type") == "CHEM_SPEC"]
        mech = next(e for e in cd if e.get("type") == "MECHANISM")
        self.third = {s["name"] for s in specs if s.get("third body")}
        self.species = [s["name"] for s in specs if not s.get("third body")]
        self.props = {s["name"]: s for s in specs}
        self.init_conc = np.array([self.props[s].get("initial concentration [mol m-3]", 0.0)
                                   for s in self.species])
        ic = d.get("initial conditions", {})
        self.T0 = float(ic.get("temperature", 298.15))
        self.P0 = float(ic.get("pressure", 101325.0))
        self.params0 = {k: float(v) for k, v in d.get("initial rate parameters", {}).items()}
        self.reactions = mech["reactions"]
        pairs = []
        for r in self.reactions:
            if r["type"] == "SURFACE":
                reactants = [r["gas-phase reactant"]]
                products = [(p["species"], float(p.get("yield", 1.0))) for p in r["gas-phase products"]]
            else:
                reactants = list(r["reactants"])
                products = [(p["species"], float(p.get("yield", 1.0))) for p in r.get("products", [])]
            pairs.append((reactants, products))
        self.mech = Mechanism(self.species, pairs, third_bodies=self.third)
        self.index = {s: i for i, s in enumerate(self.species)}
        self.photolysis_reactions = [j for j, r in enumerate(self.reactions)
                                     if r["type"] == "USER_DEFINED" and r["MUSICA name"].startswith("j")]

    # ---- rate constants -------------------------------------------------------
    def rate_constants(self, T, P, params):
        """k of shape (C, R) for temperatures T, pressures P (C,) and rate
        parameters params: name -> (C,)."""
        T = np.asarray(T, dtype=float)
        P = np.asarray(P, dtype=float)
        C = len(T)
        M = P / (R_GAS * T)
        k = np.empty((C, len(self.reactions)))
        for j, r in enumerate(self.reactions):
            t = r["type"]
            if t == "ARRHENIUS":
                A = r.get("A", 1.0); B = r.get("B", 0.0); Cc = r.get("C", 0.0)
                D = r.get("D", 300.0); E = r.get("E", 0.0)
                k[:, j] = A * np.exp(Cc / T) * (T / D) ** B * (1.0 + E * P)
            elif t == "TROE":
                k0 = r.get("k0_A", 1.0) * np.exp(r.get("k0_C", 0.0) / T) * (T / 300.0) ** r.get("k0_B", 0.0)
                kinf = r.get("kinf_A", 1.0) * np.exp(r.get("kinf_C", 0.0) / T) * (T / 300.0) ** r.get("kinf_B", 0.0)
                Fc = r.get("Fc", 0.6); N = r.get("N", 1.0)
                ratio = k0 * M / kinf
                k[:, j] = k0 * M / (1.0 + ratio) * Fc ** (N / (N + np.log10(ratio) ** 2))
            elif t == "USER_DEFINED":
                k[:, j] = params[r["MUSICA name"]] * r.get("scaling factor", 1.0)
            elif t == "SURFACE":
                label = r["MUSICA name"]
                radius = params[label + ".effective radius [m]"]
                number = params[label + ".particle number concentration [# m-3]"]
                sp = self.props[r["gas-phase reactant"]]
                Dc = sp["diffusion coefficient [m2 s-1]"]
                mw = sp["molecular weight [kg mol-1]"]
                v = np.sqrt(8.0 * R_GAS / (math.pi * mw) * T)
                gamma = r.get("reaction probability", 1.0)
                k[:, j] = 4.0 * number * math.pi * radius * radius / (radius / Dc + 4.0 / (v * gamma))
            else:
                raise ValueError(f"rate type {t} not in TS1")
        power = self.mech.third_body_power
        if np.any(power):
            k *= M[:, None] ** power[None, :]
        return k


def read_csv_table(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[float(v) for v in row] for row in reader if row]
    return header, np.array(rows)


class Cells:
    """The cells of a CELLS_CSV: TS1's initial state, then each row's values."""

    def __init__(self, ts1: TS1, path, which=None):
        header, rows = read_csv_table(path)
        if which is not None:
            rows = rows[np.asarray(which)]
        C = rows.shape[0]
        self.C = C
        self.T = np.full(C, ts1.T0)
        self.P = np.full(C, ts1.P0)
        self.Y = np.tile(ts1.init_conc, (C, 1))
        self.params = {k: np.full(C, v) for k, v in ts1.params0.items()}
        for col, name in enumerate(header):
            if name == "temperature":
                self.T = rows[:, col].copy()
            elif name == "pressure":
                self.P = rows[:, col].copy()
            elif name in ts1.index:
                self.Y[:, ts1.index[name]] = rows[:, col]
            else:
                self.params[name] = rows[:, col].copy()

    def subset(self, sl):
        other = object.__new__(Cells)
        other.C = len(range(*sl.indices(self.C)))
        other.T = self.T[sl].copy()
        other.P = self.P[sl].copy()
        other.Y = self.Y[sl].copy()
        other.params = {k: v[sl].copy() for k, v in self.params.items()}
        return other


class Scenario:
    def __init__(self, path):
        self.path = Path(path)
        values = {}
        for line in self.path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            key, value = (s.strip() for s in line.split("=", 1))
            values[key] = value
        base = self.path.parent
        self.name = base.name
        self.cells = base / values["cells"]
        self.dt = float(values["dt"])
        self.steps = int(values["steps"])
        self.snapshot_every = int(values.get("snapshot_every", 0))
        _, table = read_csv_table(base / values["photolysis"])
        self._ft = table[:, 0]
        self._fv = table[:, 1]
        self.reference_snapshots = base / "reference_snapshots.csv"
        self.reference_final = base / "reference_final.csv"

    def factor(self, t):
        return float(np.interp(t, self._ft, self._fv))

    def snapshot_steps(self):
        if self.snapshot_every <= 0:
            return [self.steps]
        steps = list(range(self.snapshot_every, self.steps + 1, self.snapshot_every))
        if not steps or steps[-1] != self.steps:
            steps.append(self.steps)
        return steps


def run_scenario(scenario: Scenario, solver, cells: Cells, ts1: TS1, rtol, atol,
                 group_size=32, constrained=(), warmup_steps=0, steps=None, on_step=None):
    """Run a solver on a scenario's cells in groups that share a step. Returns
    (snapshots: step -> (C, n), counts, statuses)."""
    n = ts1.mech.n
    atol_v = np.full(n, atol) if np.isscalar(atol) else np.asarray(atol)
    mass = np.ones(n)
    for s in constrained:
        mass[ts1.index[s]] = 0.0
    nsteps = scenario.steps if steps is None else steps
    snap_steps = [s for s in scenario.snapshot_steps() if s <= nsteps] if steps is None else [nsteps]
    snapshots = {s: np.empty((cells.C, n)) for s in snap_steps}
    counts = Counts()
    statuses = []
    photolysis = [k for k in cells.params if k.startswith("j")]
    for g0 in range(0, cells.C, group_size):
        sl = slice(g0, min(g0 + group_size, cells.C))
        group = cells.subset(sl)
        base = {k: group.params[k].copy() for k in photolysis}
        Y = group.Y.copy()
        memory = solver.new_memory()
        k_static = None
        for step in range(-warmup_steps, nsteps):
            f = scenario.factor(max(step, 0) * scenario.dt)
            for name in photolysis:
                group.params[name] = base[name] * f
            if k_static is None:
                k_static = ts1.rate_constants(group.T, group.P, group.params)
                k = k_static
            else:
                k = k_static.copy()
                for j in ts1.photolysis_reactions:
                    r = ts1.reactions[j]
                    k[:, j] = group.params[r["MUSICA name"]] * r.get("scaling factor", 1.0)
            system = ChemSystem(ts1.mech, k, mass=mass, dense=False)
            Y, status = solver.solve(system, Y, scenario.dt, atol_v, rtol, counts, memory)
            if status != "converged":
                statuses.append((g0, step, status))
            if on_step is not None:
                on_step(g0, step, Y)
            if step + 1 in snapshots:
                snapshots[step + 1][sl] = Y
    return snapshots, counts, statuses


def score(result: dict, cells_index, species, reference_path, threshold=1e-14, steps=None):
    """Digits overall and per species (Sandu et al.), as the scenarios' score.py
    computes them, over the given cells: result maps step -> (C, n) in the order
    of cells_index. Returns (overall, per-species dict, gate)."""
    header, rows = read_csv_table(reference_path)
    keyed = header[:2] == ["step", "cell"]
    start = 4 if keyed else 2
    ref_species = header[start:]
    col = {s: start + i for i, s in enumerate(ref_species)}
    sp_index = {s: i for i, s in enumerate(species)}
    want = list(cells_index)
    pos = {c: i for i, c in enumerate(want)}
    sq_by_species = {s: [] for s in ref_species}
    gate = True
    for step, Y in result.items():
        if not np.all(np.isfinite(Y)) or np.any(Y < -1e-20):
            gate = False
        if keyed:
            sel = rows[(rows[:, 0] == step) & np.isin(rows[:, 1], want)]
            pairs = [(Y[pos[int(row[1])]], row) for row in sel]
        else:
            pairs = [(Y[i], rows[c]) for i, c in enumerate(want)]
        for y_cell, row in pairs:
            for s in ref_species:
                ref = row[col[s]]
                if ref < threshold or s not in sp_index:
                    continue
                y = y_cell[sp_index[s]]
                sq_by_species[s].append(((y - ref) / ref) ** 2)
    all_sq = [v for vs in sq_by_species.values() for v in vs]

    def digits(values):
        if not values:
            return float("nan")
        ms = float(np.mean(values))
        return 16.0 if ms <= 0.0 else min(16.0, -math.log10(math.sqrt(ms)))

    per = {s: digits(v) for s, v in sq_by_species.items() if v}
    return digits(all_sq), per, gate
