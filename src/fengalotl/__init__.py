from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("fengalotl")
except PackageNotFoundError:
    __version__ = "0.0.0_fallback"
