"""Export this nonmagnetic NbO SCF snapshot using only Python's standard library."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parent
SOURCE = None
OUTPUT = BASE / "data"
CANDIDATES = (25, 26, 27)  # VASP's one-based, energy-ordered indices


def write_csv(name, rows, fields):
    with (OUTPUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    global SOURCE, OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="Directory containing original VASP SCF outputs (not bundled).")
    parser.add_argument("--output", type=Path, default=BASE / "data")
    args = parser.parse_args()
    SOURCE, OUTPUT = args.source.resolve(), args.output.resolve()
    root = ET.parse(SOURCE / "vasprun.xml").getroot()
    parameters = root.find("parameters")
    def parameter(name):
        return parameters.find(f".//i[@name='{name}']").text.strip()

    assert int(parameter("ISPIN")) == 1, "This exporter assumes ISPIN=1."
    for name in ("LSORBIT", "LNONCOLLINEAR"):
        value = parameters.find(f".//i[@name='{name}']")
        assert value is None or value.text.strip().upper() in ("F", ".FALSE."), name
    calc = root.findall("calculation")[-1]
    ef = float(calc.find("dos/i[@name='efermi']").text)
    kpoints = [[float(x) for x in v.text.split()]
               for v in root.findall("kpoints/varray[@name='kpointlist']/v")]
    weights = [float(v.text) for v in root.findall("kpoints/varray[@name='weights']/v")]
    # For this ISPIN=1 XML, occupations are per spin (0..1).
    eigen = [[[float(x) for x in row.text.split()] for row in block.findall("r")]
             for block in calc.find("eigenvalues/array/set/set").findall("set")]
    nk, nb = len(eigen), len(eigen[0])
    assert len(kpoints) == len(weights) == nk
    assert all(len(block) == nb for block in eigen)
    assert abs(sum(weights) - 1) < 1e-5
    assert all(0 <= row[1] <= 1 for block in eigen for row in block)
    nelect = float(parameter("NELECT"))
    integrated_nelect = 2 * sum(weights[k] * sum(row[1] for row in eigen[k])
                                 for k in range(nk))
    assert abs(integrated_nelect - nelect) < 1e-3

    # Independent check against EIGENVAL (XML energies have fewer printed digits).
    lines = (SOURCE / "EIGENVAL").read_text().splitlines()
    eig_ne, eig_nk, eig_nb = map(int, lines[5].split())
    assert (eig_ne, eig_nk, eig_nb) == (int(nelect), nk, nb)
    cursor, max_diff = 6, 0.0
    for k in range(nk):
        while not lines[cursor].strip():
            cursor += 1
        coords_weight = list(map(float, lines[cursor].split()))
        assert max(abs(a - b) for a, b in zip(coords_weight[:3], kpoints[k])) < 1e-7
        assert abs(coords_weight[3] - weights[k]) < 1e-7
        cursor += 1
        for b in range(nb):
            fields = lines[cursor].split()
            cursor += 1
            assert int(fields[0]) == b + 1
            max_diff = max(max_diff, abs(float(fields[1]) - eigen[k][b][0]))
    assert max_diff < 1e-4

    # Sum PAW projections by species and angular momentum, without interpreting
    # them as full-wavefunction probabilities. POSCAR's order must be Nb3 O3.
    poscar = (SOURCE / "POSCAR").read_text().splitlines()
    assert poscar[5].split() == ["Nb", "O"]
    assert list(map(int, poscar[6].split())) == [3, 3]
    projections = {}
    ik = ib = None
    with (SOURCE / "PROCAR").open() as stream:
        for line in stream:
            match = re.match(r"\s*k-point\s+(\d+)\s*:", line)
            if match:
                ik = int(match.group(1)) - 1
            match = re.match(r"\s*band\s+(\d+)\s*#", line)
            if match:
                ib = int(match.group(1)) - 1
            if line.split()[:1] != ["ion"]:
                continue
            labels = line.split()[1:]
            totals = {f"paw_{species}_{orb}": 0.0
                      for species in ("Nb", "O") for orb in ("s", "p", "d")}
            for atom in range(1, 7):
                row = next(stream).split()
                assert int(row[0]) == atom
                species = "Nb" if atom <= 3 else "O"
                for label, value in zip(labels, map(float, row[1:])):
                    if label == "tot":
                        continue
                    group = "d" if label.startswith("d") or label == "x2-y2" else label[0]
                    assert group in ("s", "p", "d")
                    totals[f"paw_{species}_{group}"] += value
            assert (ik, ib) not in projections
            projections[ik, ib] = totals
    assert len(projections) == nk * nb

    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows, summary = [], []
    for k, (coords, weight) in enumerate(zip(kpoints, weights)):
        for b, (energy, occupation) in enumerate(eigen[k]):
            row = dict(k_index=k + 1, band_index=b + 1,
                       k1=coords[0], k2=coords[1], k3=coords[2],
                       weight=weight, energy_eV=energy, xi_eV=energy - ef,
                       occupation_per_spin=occupation)
            row.update(projections[k, b])
            row["paw_sum"] = sum(projections[k, b].values())
            rows.append(row)
    fields = list(rows[0])
    candidates = [row for row in rows if row["band_index"] in CANDIDATES]
    write_csv("candidate_bands_25_27.csv", candidates, fields)
    for b in range(1, nb + 1):
        band = [row for row in rows if row["band_index"] == b]
        lo, hi = min(row["xi_eV"] for row in band), max(row["xi_eV"] for row in band)
        closest = min(band, key=lambda row: abs(row["xi_eV"]))
        item = dict(band_index=b, min_xi_eV=lo, max_xi_eV=hi,
                    crosses_EF_sampled=lo < 0 < hi,
                    electrons_per_cell=2 * sum(row["weight"] * row["occupation_per_spin"]
                                               for row in band),
                    closest_k_index=closest["k_index"], closest_xi_eV=closest["xi_eV"])
        for window in (0.05, 0.1, 0.2, 0.5):
            item[f"n_irreducible_states_abs_xi_lt_{window}_eV"] = sum(
                abs(row["xi_eV"]) < window for row in band)
        summary.append(item)
    write_csv("band_summary.csv", summary, list(summary[0]))
    metadata = dict(
        source_description="User-computed NbO PBE static SCF; original VASP outputs archived separately",
        material_reference="https://legacy.materialsproject.org/materials/mp-2311/",
        efermi_eV=ef, nelect=nelect, integrated_nelect=integrated_nelect,
        nk_irreducible=nk, nbands=nb, weight_sum=sum(weights),
        ispin=1, soc=False, mesh=[16, 16, 16],
        candidate_bands=list(CANDIDATES),
        sampled_crossing_bands=[s["band_index"] for s in summary if s["crosses_EF_sampled"]],
        all_band_rows_checked=len(rows), candidate_rows=len(candidates),
        k_coordinate_convention="fractional reciprocal-lattice coordinates",
        energy_reference="energy_eV is VASP eigenvalue; xi_eV = energy_eV - efermi_eV",
        occupation_convention="XML ISPIN=1 occupation per spin; multiply by 2 for electron count",
        projection_convention="raw PAW species/angular-momentum sums; not full probabilities",
        eigenval_xml_max_energy_difference_eV=max_diff,
        interaction_available=False,
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(SOURCE.iterdir()) if p.is_file()},
    )
    (OUTPUT / "metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Checked {len(rows)} states; exported {len(candidates)} candidate states.")
    print(f"EF={ef:.8f} eV; electrons/cell={integrated_nelect:.8f}")
    print("Crossing bands:", metadata["sampled_crossing_bands"])
    print(f"EIGENVAL/XML maximum difference={max_diff:.8g} eV")
    print("Output:", OUTPUT)


if __name__ == "__main__":
    main()
