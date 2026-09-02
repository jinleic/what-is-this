"""Boot stub: register the shared d4core dependency under the name that
Campaign A's gate_corr_env.py expects in sys.modules BEFORE exec_module,
exactly as A/B/C did (their __file__-relative path resolves to
campaigns/20260831T082425Z_kg_direct_d4_restart/code/core.py; that byte
target exists unchanged, and its import here re-creates their module
construction pattern for this run's byte-shared instrument).
"""
import importlib.util as ilu
import os
import sys

sys.dont_write_bytecode = True
_CODE = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.join(_CODE, "..", "..", "20260831T082425Z_kg_direct_d4_restart", "code", "core.py")

_spec = ilu.spec_from_file_location("d4core", CORE)
_d4core = ilu.module_from_spec(_spec)
sys.modules["d4core"] = _d4core
_spec.loader.exec_module(_d4core)
