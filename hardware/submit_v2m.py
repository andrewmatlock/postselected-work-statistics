"""
submit_v2m.py — registered submission wrapper for PREREGISTRATION_witness_v2m.

Waits for ibm_marrakesh to leave maintenance, re-evaluates the registered
health rule on the then-current published calibration, and only then runs
hardware/run_witness.py (unmodified). Polls every 10 min; after a failed
health check, re-checks hourly (calibration updates ~daily). Gives up after
72 h with a clear message.
"""
import subprocess, sys, time

from qiskit_ibm_runtime import QiskitRuntimeService

DEADLINE = time.time() + 72 * 3600
service = QiskitRuntimeService()

def health(b):
    p = b.properties(refresh=True)
    ro = [p.readout_error(q) for q in (53, 54, 55)]
    cz = [b.target["cz"][(53, 54)].error, b.target["cz"][(54, 55)].error]
    ok = all(e <= 0.03 for e in ro) and all(e <= 0.01 for e in cz)
    print(f"health: readout {['%.4f' % e for e in ro]} cz {['%.4f' % e for e in cz]}"
          f" -> {'PASS' if ok else 'FAIL'}", flush=True)
    return ok

while time.time() < DEADLINE:
    try:
        b = service.backend("ibm_marrakesh")
        st = b.status()
        print(f"{time.strftime('%H:%M:%S')} status={st.status_msg} pending={st.pending_jobs}",
              flush=True)
        if st.status_msg != "maintenance" and st.operational:
            if health(b):
                print("submitting...", flush=True)
                r = subprocess.run([sys.executable, "-u", "run_witness.py",
                                    "hardware", "ibm_marrakesh"])
                sys.exit(r.returncode)
            time.sleep(3600)
            continue
    except Exception as e:
        print(f"transient: {type(e).__name__}: {e}", flush=True)
    time.sleep(600)

print("gave up after 72 h: backend never passed status+health; resubmit later",
      flush=True)
sys.exit(2)
