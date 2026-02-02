"""
moducomp: metabolic module completeness and complementarity for microbiomes.
"""

__version__ = "0.7.9"
__author__ = "Juan C. Villada"
__email__ = "jvillada@lbl.gov"
__title__ = "moducomp"
__description__ = "Metabolic module completeness and complementarity for microbiomes"
__license__ = "BSD-3-Clause"

from .moducomp import app

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__title__",
    "__description__",
    "__license__",
    "app",
]


def get_version() -> str:
    """Return the package version."""
    return __version__


def main() -> None:
    """Entry point for the moducomp command line interface."""
    app()
