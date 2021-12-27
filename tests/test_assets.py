"""Test that the assest we distribute are actually present."""
from subprocess import run
import gzip
import os
from pyprojroot import here
from seqlike.version import __version__

os.chdir(here())
print(os.getcwd())


def test_free_mono_font_exists():
    """Tests that the FreeMono.ttf font exists in source distribution."""
    run("rm -r dist/*".split(" "), check=False)
    run("python -m build".split(" "), check=True)
    dist_dir = here() / "dist"
    os.chdir(dist_dir)
    run(f"tar -xvf seqlike-{__version__}.tar.gz".split(" "), check=True)
    file_path = dist_dir / f"seqlike-{__version__}/seqlike/FreeMono.ttf"

    assert file_path.exists()

    run("rm -rf dist/*".split(" "), check=False)
