#!/usr/bin/env python3
"""Sanity checks for calibration-hit simulation output.

  python3 check_calibhits.py calib.edm4hep.root [--baseline baseline.edm4hep.root]

Per calorimeter and event it reports the active-layer energy, the calibration
(total deposited) energy, and their ratio. It also checks that:
  * every active cell also has a calibration hit;
  * with --baseline (same seed, steer_baseline.py), the active collections are
    hit-for-hit identical to the standard simulation.
Run inside the mucoll container (needs podio + edm4hep python bindings).
"""
import argparse
import sys

from podio.root_io import Reader

CALOS = ["ECalBarrelCollection", "ECalEndcapCollection",
         "HCalBarrelCollection", "HCalEndcapCollection"]


def hit_map(coll):
    return {h.getCellID(): h for h in coll}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("calib", help="ddsim output produced with steer_calib.py")
    ap.add_argument("--baseline", help="same-seed output produced with steer_baseline.py")
    ap.add_argument("--suffix", default="Calib")
    ap.add_argument("--max-events", type=int, default=-1)
    args = ap.parse_args()

    events = Reader(args.calib).get("events")
    base_events = Reader(args.baseline).get("events") if args.baseline else None
    n_fail = 0
    totals = {c: [0.0, 0.0] for c in CALOS}

    for iev, ev in enumerate(events):
        if 0 <= args.max_events <= iev:
            break
        print(f"--- event {iev}")
        bev = base_events[iev] if base_events is not None else None
        for name in CALOS:
            active = hit_map(ev.get(name))
            calib = hit_map(ev.get(name + args.suffix))
            e_act = sum(h.getEnergy() for h in active.values())
            e_cal = sum(h.getEnergy() for h in calib.values())
            totals[name][0] += e_act
            totals[name][1] += e_cal

            missing = [c for c in active if c not in calib]
            if missing:
                n_fail += 1
                print(f"  FAIL {name}: {len(missing)} active cells without a calib hit")
            n_contrib_cal = sum(h.getContributions().size() for h in calib.values())
            only_passive = sum(h.getEnergy() for c, h in calib.items() if c not in active)
            ratio = e_cal / e_act if e_act > 0 else float("nan")
            print(f"  {name:22s} active: {len(active):6d} hits {e_act:9.4f} GeV | "
                  f"calib: {len(calib):6d} hits {e_cal:9.4f} GeV ({n_contrib_cal} contribs) | "
                  f"calib/active = {ratio:7.2f} | E in passive-only cells = {only_passive:.4f} GeV")

            if bev is not None:
                ref = hit_map(bev.get(name))
                bad = 0
                if set(ref) != set(active):
                    bad += 1
                    print(f"  FAIL {name}: cellID sets differ from baseline "
                          f"({len(set(ref) ^ set(active))} cells)")
                for c, h in ref.items():
                    a = active.get(c)
                    if a is None:
                        continue
                    if (abs(a.getEnergy() - h.getEnergy()) > 1e-9 * max(1.0, h.getEnergy())
                            or a.getContributions().size() != h.getContributions().size()):
                        bad += 1
                if bad:
                    n_fail += 1
                    print(f"  FAIL {name}: {bad} differences w.r.t. baseline")
                else:
                    print(f"  OK   {name}: identical to baseline ({len(ref)} hits)")

    print("--- totals over events")
    for name, (ea, ec) in totals.items():
        r = ec / ea if ea > 0 else float("nan")
        print(f"  {name:22s} active {ea:10.4f} GeV  calib {ec:10.4f} GeV  ratio {r:7.2f}")
    print("RESULT:", "FAIL" if n_fail else "OK")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
