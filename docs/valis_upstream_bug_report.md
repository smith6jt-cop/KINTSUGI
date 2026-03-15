# Bug Report: `UnboundLocalError` in `registration.py:register_micro()` (v1.2.0)

**Repository:** https://github.com/MathOnco/valis
**Version:** valis-wsi 1.2.0
**Severity:** Crash — `UnboundLocalError` when calling `register_micro()` with a class and params

## Description

`Valis.register_micro()` in `registration.py` (lines 4890-4897) has a bug where passing an uninstantiated class with parameters causes an `UnboundLocalError`.

## Code with Bug

```python
# registration.py lines 4890-4897
if non_rigid_registrar_cls is not None:
    if isinstance(non_rigid_registrar_cls, type):
        if non_rigid_reg_params is None:
            non_rigid_registrar_obj = non_rigid_registrar_cls()
        else:
            nr = non_rigid_registrar_obj(**non_rigid_reg_params)  # BUG
    else:
        non_rigid_registrar_obj = non_rigid_registrar_cls
```

Two issues on line 4896:
1. `non_rigid_registrar_obj` is **not yet defined** at this point — it's only assigned in the `if non_rigid_reg_params is None` branch above. This causes `UnboundLocalError`.
2. The result is assigned to `nr` instead of `non_rigid_registrar_obj`, so even if it didn't crash, the variable used on subsequent lines (`non_rigid_registrar_obj`) would still be undefined.

## Expected Fix

```python
if non_rigid_registrar_cls is not None:
    if isinstance(non_rigid_registrar_cls, type):
        if non_rigid_reg_params is None:
            non_rigid_registrar_obj = non_rigid_registrar_cls()
        else:
            non_rigid_registrar_obj = non_rigid_registrar_cls(**non_rigid_reg_params)
    else:
        non_rigid_registrar_obj = non_rigid_registrar_cls
```

## Reproduction

```python
from valis import registration, non_rigid_registrars

registrar = registration.Valis(
    src_dir="/path/to/edf",
    dst_dir="/path/to/output",
    non_rigid_registrar_cls=None,  # Disable coarse NR
)

registrar.register()

# This call triggers the bug:
registrar.register_micro(
    max_non_rigid_registration_dim_px=4096,
    non_rigid_registrar_cls=non_rigid_registrars.OpticalFlowWarper,  # a CLASS, not instance
    non_rigid_reg_params={"smoothing_method": "gauss", "sigma_ratio": 0.01},
)
# -> UnboundLocalError: local variable 'non_rigid_registrar_obj' referenced before assignment
```

## Note

The same pattern is correctly implemented in `serial_non_rigid.py:register()` (lines 995-1001), where it properly handles both class and instance inputs. The `register_micro()` method appears to have been partially updated.

## Context

We use `register_micro()` for multi-cycle immunofluorescence registration (CODEX/mIF), passing `OpticalFlowWarper` as a class with smoothing parameters. This is the standard usage pattern documented in VALIS tutorials.

## Environment

- valis-wsi 1.2.0
- Python 3.10-3.12
- Linux (Ubuntu 24.04 / HiPerGator HPC)
