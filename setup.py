import os
import requests
from setuptools import setup, find_packages
from setuptools.command.install import install
from distutils.command.install import install as orig_install
import subprocess


class InstallWrapper(install):
    """Install MAFFT binary executables and other requirements.

    The order of priority in which we try:

    - Try to install via dpkg.
    - Try to install via brew.
    - Raise warning to install on their own.

    :sa: https://www.anomaly.net.au/blog/running-pre-and-post-install-jobs-for-your-python-packages/
    :sa: https://stackoverflow.com/questions/21915469/python-setuptools-install-requires-is-ignored-when-overriding-cmdclass
    """

    def run(self):
        # install the python requirements from `install_requires`
        # problem with setuptools.command.install and wheel setup
        # :sa: https://github.com/automl/random_forest_run/commit/f975f947979cd1bc46eca30833f8267c2f71514c
        # install.do_egg_install(self)
        orig_install.run(self)
        # install the non-python binaries by Debian
        result = dpkg_install()
        if result.returncode != 0:
            result = brew_install()
        if result.returncode != 0:
            warning_install()


def warning_install():
    """Raise a warning that none of the previous installation methods worked."""
    import warnings

    warnings.warn(
        "Installation by dpkg and brew failed. "
        "Please manually install MAFFT to enjoy SeqLike's Series Accessor's "
        "method-chained alignment capabilities. \n\n"
        "MAFFT installation instructions can be found on "
        "https://mafft.cbrc.jp/alignment/software/."
    )


def brew_install():
    """Attempt installation by brew."""
    result = subprocess.run(["brew", "update"])
    if result.returncode == 0:
        result = subprocess.run(["brew", "install", "mafft"])
    return result


def dpkg_install():
    """Install via dpkg."""
    result = subprocess.run(["apt-get", "update"])
    if result.returncode == 0:
        result = subprocess.run(["apt-get", "install", "-y", "ghostscript"])
    if result.returncode == 0:
        result = download_and_install_mafft(kind="deb")
    return result


def yum_install():
    """Placeholder for installation via yum.

    NOTE: Though I've used RedHat and CentOS, I'm not familiar enough with the commands.
    (@ericmjl)
    """
    pass


def download_and_install_mafft(kind="deb"):
    """Install MAFFT as a Debian package.

    Default to Debian install.
    """
    mafft_filenames = dict(rpm="mafft-7.487-gcc_fc6.x86_64.rpm", deb="mafft_7.487-1_amd64.deb")
    mafft_file = mafft_filenames[kind]
    deb_filepath = "/tmp/" + mafft_file

    # download the package from mafft.cbrc.jp
    if not os.path.isfile(deb_filepath):
        url = "https://mafft.cbrc.jp/alignment/software/" + mafft_file
        print("Downloading %s..." % url)
        r = requests.get(url, allow_redirects=True)
        with open(deb_filepath, "wb") as f:
            f.write(r.content)

    # install the MAFFT package by dpkg
    print("Installing %s..." % deb_filepath)
    return subprocess.run(["dpkg", "--install", deb_filepath])


setup(
    name="seqlike",
    version="1.1.2",
    packages=find_packages(),  # https://stackoverflow.com/a/22442340
    cmdclass={"install": InstallWrapper},
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "scikit-learn",
        "weblogo",
        "Pillow",
        # :note: pytest should be under tests_require, but this doesn't seem to work
        "pytest-regtest",
        "pytest",
        "multipledispatch",
    ],
    # tests_require=[
    #    'pytest-regtest',
    #    'pytest',
    # ],
)
