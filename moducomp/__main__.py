"""
Entry point for running moducomp as a module.

This allows running moducomp with:
python -m moducomp [commands]
"""

from .moducomp import app

if __name__ == "__main__":
    app()
