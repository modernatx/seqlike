import os
import subprocess
import warnings
from distutils.command.install import install as orig_install
from pathlib import Path

import requests
from setuptools import find_packages, setup
from setuptools.command.install import install

SUCCESS_RETURN_CODE = 0


class InstallWrapper(install):
    """Install MAFFT binary executables and other requirements.

    The order of priority in which we try:

    - Determine whether no installation is necessary (no_install).
    - Try to install via dpkg (debian_install).
    - Try to install via brew (brew_install).
    - Raise warning to install on their own (warning_install).

    :sa: https://www.anomaly.net.au/blog/running-pre-and-post-install-jobs-for-your-python-packages/
    :sa: https://stackoverflow.com/questions/21915469/python-setuptools-install-requires-is-ignored-when-overriding-cmdclass
    """

    def run(self):
        # install the python requirements from `install_requires`
        # problem with setuptools.command.install and wheel setup
        # :sa: https://github.com/automl/random_forest_run/commit/f975f947979cd1bc46eca30833f8267c2f71514c
        # install.do_egg_install(self)
        orig_install.run(self)
        result = no_install()
        if result.returncode != SUCCESS_RETURN_CODE:
            result = debian_install()
        if result.returncode != SUCCESS_RETURN_CODE:
            print("Debian-based installation failed. Trying with brew instead...")
            result = brew_install()
        if result.returncode != SUCCESS_RETURN_CODE:
            warning_install()


def no_install():
    """Run `which mafft` to see if mafft is already installed."""
    result = subprocess.run(["which", "mafft"])
    return result


def brew_install():
    """Attempt installation by brew."""
    result = subprocess.run(["which", "brew"])
    if result.returncode == SUCCESS_RETURN_CODE:
        result = subprocess.run(["brew", "update"])
    if result.returncode == SUCCESS_RETURN_CODE:
        result = subprocess.run(["brew", "install", "mafft"])
    return result


def debian_install():
    """Install MAFFT on debian."""
    result = subprocess.run(["which", "apt-get"])
    if result.returncode == SUCCESS_RETURN_CODE:
        result = subprocess.run(["apt-get", "update"])
    if result.returncode == SUCCESS_RETURN_CODE:
        result = subprocess.run(["apt-get", "install", "-y", "ghostscript"])
    if result.returncode == SUCCESS_RETURN_CODE:
        deb_filepath = download_mafft(kind="deb")
        print(f"Installing {deb_filepath}...")
        result = subprocess.run(["dpkg", "--install", deb_filepath])
    return result


def rpm_install():
    """Install MAFFT using RPM."""
    raise NotImplementedError("Installation by RPM not yet supported.")


def warning_install():
    """Raise a warning that none of the previous installation methods worked."""
    warnings.warn(
        "Installation by dpkg and brew failed. "
        "Please manually install MAFFT to enjoy SeqLike's Series Accessor's "
        "method-chained alignment capabilities. \n\n"
        "MAFFT installation instructions can be found on "
        "https://mafft.cbrc.jp/alignment/software/."
    )


def download_mafft(kind="deb") -> Path:
    """Download MAFFT from the original website."""
    acceptable_kinds = ["deb", "rpm"]
    if kind not in acceptable_kinds:
        raise NameError(f"`kind` should be one of {acceptable_kinds}.")

    if kind == "deb":
        mafft_file = "mafft_7.487-1_amd64.deb"
    if kind == "rpm":
        mafft_file = "mafft-7.487-gcc_fc6.x86_64.rpm"

    deb_filepath = Path("/tmp/") / mafft_file
    if not os.path.isfile(deb_filepath):
        url = "https://mafft.cbrc.jp/alignment/software/" + mafft_file
        print("Downloading %s..." % url)
        r = requests.get(url, allow_redirects=True)
        with open(deb_filepath, "wb") as f:
            f.write(r.content)
    return deb_filepath


README = Path("README.md").read_text()

install_requires = [
    "Pillow",
    "biopython",
    "lazy-loader",
    "multipledispatch",
    "numpy",
    "pandas",
    "python-codon-tables",
    "scikit-learn",
    "weblogo",
]
test_requires = [
    "pytest",
    "pytest-regtest",
]

setup(
    name="seqlike",
    version="1.3.4",
    packages=find_packages(),  # https://stackoverflow.com/a/22442340
    cmdclass={"install": InstallWrapper},
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/modernatx/seqlike",
    project_urls={
        "Documentation": "https://modernatx.github.io/seqlike",
        "Source Code": "https://github.com/modernatx/seqlike",
        "Issue Tracker": "https://github.com/modernatx/seqlike/issues",
    },
    install_requires=install_requires,
    extra_requires={
        "notebook": ["bokeh>=3.0.2"],  # Used in ipython/jupyter
        "test": test_requires,
    },
    package_data={"seqlike": ["*.ttf"]},
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
    ],
)
