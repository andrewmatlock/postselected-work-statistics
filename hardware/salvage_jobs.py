"""
salvage_jobs.py — retrieve results from already-completed IBM jobs by ID.

Retrieving a finished job costs NO QPU time. This reassembles whatever
completed before the usage cap into the standard counts JSON.

Usage:  python salvage_jobs.py
Output: ctc_shallow_salvaged.json
"""
import json
import numpy as np
from qiskit_ibm_runtime import QiskitRuntimeService

# job IDs printed during the run, grouped by direction and input state
JOBS = {
    "forward": {
        0: ["d9mv6nuij12s73ft9efg", "d9mv6s460llc73c9g17g", "d9mv6v0qs0bc73e2qpsg",
            "d9mv73uij12s73ft9es0", "d9mv774sfqic73aqe9pg"],
        1: ["d9mv7a460llc73c9g1mg", "d9mv7d4sfqic73aqea0g", "d9mv7gc60llc73c9g1ug",
            "d9mv7jc60llc73c9g220", "d9mv7mcsfqic73aqeab0"],
        2: ["d9mv7p4sfqic73aqeaeg", "d9mv7rssfqic73aqeahg", "d9mv7v460llc73c9g2eg",
            "d9mv82460llc73c9g2hg", "d9mv86csfqic73aqear0"],
        3: ["d9mv89csfqic73aqeaug", "d9mv8c460llc73c9g2s0", "d9mv8f8qs0bc73e2qrcg",
            "d9mv8ic60llc73c9g31g", "d9mv9bc60llc73c9g3q0"],
    },
    "reverse": {
        0: ["d9mv9eeij12s73ft9h4g", "d9mv9h0qs0bc73e2qsg0", "d9mv9keij12s73ft9hb0",
            "d9mv9o0qs0bc73e2qso0", "d9mv9r4sfqic73aqech0"],
        1: ["d9mv9u4sfqic73aqeck0", "d9mva16ij12s73ft9hng", "d9mva3s60llc73c9g4m0",
            "d9mva6uij12s73ft9hug", "d9mva9ssfqic73aqecv0"],
        2: ["d9mvad460llc73c9g4vg", "d9mvag6ij12s73ft9i7g", "d9mvaioqs0bc73e2qtlg"],
        # input 3 never ran — cap hit
    },
}

service = QiskitRuntimeService()
out = {"forward": {}, "reverse": {}}
totals = {"forward": {}, "reverse": {}}

for direction, groups in JOBS.items():
    for i, ids in groups.items():
        acc, got, shots = {}, 0, 0
        for jid in ids:
            try:
                job = service.job(jid)
                st = job.status()
                st = st if isinstance(st, str) else st.name
                if st not in ("DONE", "COMPLETED"):
                    print(f"  {direction[0].upper()}{i} {jid}: {st} — skipped")
                    continue
                counts = job.result()[0].data.c.get_counts()
                for k, v in counts.items():
                    acc[k] = acc.get(k, 0) + int(v)
                got += 1; shots += sum(counts.values())
            except Exception as e:
                print(f"  {direction[0].upper()}{i} {jid}: retrieval failed ({type(e).__name__})")
        if acc:
            out[direction][str(i)] = acc
            totals[direction][i] = shots
            print(f"{direction[0].upper()}{i}: {got}/{len(ids)} jobs, {shots} shots")
        else:
            print(f"{direction[0].upper()}{i}: NO DATA")

json.dump({"forward": out["forward"], "reverse": out["reverse"],
           "backend": "ibm_marrakesh", "shots_per_input": min(
               [v for d in totals.values() for v in d.values()] or [0]),
           "circuit": "shallow_4cx", "partial": True, "shot_totals": totals},
          open("ctc_shallow_salvaged.json", "w"), indent=2)
print("\nwrote ctc_shallow_salvaged.json")
print("inputs recovered — forward:", sorted(out["forward"].keys()),
      " reverse:", sorted(out["reverse"].keys()))
