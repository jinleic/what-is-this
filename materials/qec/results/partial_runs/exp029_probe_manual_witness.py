import json
import sys
import stim
sys.path.insert(0, 'experiments')
import exp029_circuit_distance as e

c, d, b = e.build_exact_circuit(sys.argv[1] if len(sys.argv) > 1 else 'gross')
if b['key'] == 'gross':
    cert = json.load(open('results/certificates/bb_distance_144_12_12.json'))
    support = cert['witness_X']['support']
else:
    support = [7, 21, 35, 43, 57, 71, 77, 85, 99, 113, 121, 135]
explained = c.explain_detector_error_model_errors(reduce_to_one_representative_error=False)
by_q = {}
for er in explained:
    dm = om = 0
    for t in er.dem_error_terms:
        x = t.dem_target
        if x.is_relative_detector_id(): dm ^= 1 << x.val
        elif x.is_logical_observable_id(): om ^= 1 << x.val
    for loc in er.circuit_error_locations:
        text = str(loc)
        if '(after 0 TICKs)' not in text or 'X_ERROR(0.001)' not in text:
            continue
        pp = loc.flipped_pauli_product
        if len(pp) != 1:
            continue
        target = pp[0].gate_target
        if not target.is_x_target:
            continue
        q = target.value
        by_q[q] = (dm, om, er)
print('mapped qubits', sorted(by_q))
print('mapped', len(by_q), 'support', support)
dx = ox = 0
for q in support:
    dm, om, er = by_q[q]
    dx ^= dm
    ox ^= om
    print(q, dm.bit_count(), e.mask_ids(om))
print('xor detectors', dx.bit_count(), e.mask_ids(dx))
print('xor obs', ox.bit_count(), e.mask_ids(ox))
