from __future__ import annotations

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
FROZEN_ROOT = PACKAGE_ROOT

LAYER_FILES = {
    "CCAM": "Changes in color of aerial mycelium (CCAM).xlsx",
    "CCVM": "Changes in color of vegetative mycelium (CCVM).xlsx",
    "CMP": "Changes in mycelium production (CMP).xlsx",
    "CS": "Changes in sporulation (CS).xlsx",
    "IG": "Inhibition of growth (IG).xlsx",
    "IC": "Invasion of colony (IC).xlsx",
    "RP": "Release of a pigment (RP).xlsx",
    "IRP": "Inhibition of release of pigment (IRP).xlsx",
    "IAC_RAC": "Induction of antimicrobial compounds of the other (IAC) release of antibacterial compounds (RAC).xlsx",
    "IAC_RDE": "Induction of antimicrobial compounds of the other (IAC) release of degrading enzyme (RDE).xlsx",
    "IRAC": "Inhibition of release of antimicrobial compounds (IRAC).xlsx",
    "RAC": "Release of antimicrobial compounds (RAC).xlsx",
}


def find_layer_file(filename: str) -> Path:
    path = PACKAGE_ROOT / 'data' / 'matrices_xls' / filename
    if path.exists():
        return path
    raise FileNotFoundError(filename)


def frozen_notebook_path() -> Path:
    return PACKAGE_ROOT / 'root' / '_notebook.py'
