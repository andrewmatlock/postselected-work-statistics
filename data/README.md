# Data

Raw measurement counts from quantum hardware runs, as JSON:
`{"forward": {input: {bitstring: count}}, "reverse": {...}, "backend": str, "shots_per_input": int}`

Bitstrings are Qiskit-ordered: leftmost character is the ancilla (herald) bit.
Place `ctc_hardware_counts.json`, `sim_counts.json`, and `ctc_shallow_salvaged.json` here.
