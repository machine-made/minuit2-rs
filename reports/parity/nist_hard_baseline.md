# NIST hard-dataset baseline (iminuit vs minuit2-rs)

Priority genuine-gap targets: none

ROOT runner: skipped; upstream iminuit is the C++ Minuit2 comparator for this bead.

| Dataset | iminuit s1 | iminuit s2 | minuit2-rs s1 | minuit2-rs s2 |
| --- | --- | --- | --- | --- |
| Lanczos3 | FAIL (valid=True, fval=3.3808e-05, params=False, nfcn=337) | FAIL (valid=True, fval=5.1388e-05, params=False, nfcn=226) | FAIL (valid=True, fval=4.25361e-05, params=False, nfcn=424) | FAIL (valid=True, fval=4.25209e-05, params=False, nfcn=442) |
| BoxBOD | OK (valid=True, fval=1168.01, params=True, nfcn=68) | OK (valid=True, fval=1168.01, params=True, nfcn=82) | OK (valid=True, fval=1168.01, params=True, nfcn=73) | OK (valid=True, fval=1168.01, params=True, nfcn=73) |
| MGH09 | FAIL (valid=True, fval=0.000494688, params=False, nfcn=177) | FAIL (valid=True, fval=0.000721538, params=False, nfcn=158) | FAIL (valid=True, fval=0.000341093, params=False, nfcn=275) | FAIL (valid=True, fval=0.000340651, params=False, nfcn=293) |
| Hahn1 | OK (valid=True, fval=1.53244, params=True, nfcn=581) | FAIL (valid=True, fval=44997.1, params=False, nfcn=580) | OK (valid=True, fval=1.53245, params=True, nfcn=474) | FAIL (valid=True, fval=36.4401, params=False, nfcn=1478) |

- **Lanczos3**: parity-failure. Certified residual SS=1.61172e-08; parameter tolerance=1e-02. iminuit does not reach the certified solution from Start 2; minuit2-rs does not reach it under the same strategies.
- **BoxBOD**: parity-success. Certified residual SS=1168.01; parameter tolerance=1e-02. iminuit reaches the certified solution from Start 2; minuit2-rs reaches it under the same strategies.
- **MGH09**: parity-failure. Certified residual SS=0.000307506; parameter tolerance=1e-02. iminuit does not reach the certified solution from Start 2; minuit2-rs does not reach it under the same strategies.
- **Hahn1**: parity-success. Certified residual SS=1.53244; parameter tolerance=1e-02. iminuit reaches the certified solution from Start 2; minuit2-rs reaches it under the same strategies.
