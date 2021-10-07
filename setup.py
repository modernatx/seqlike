import os
import requests
from setuptools import setup, find_packages
from setuptools.command.install import install
from distutils.command.install import install as orig_install
import subprocess


class InstallWrapper(install):
    """Install MAFFT binary executables and other requirements.
    Assumes Debian linux; has not been tested with other linux flavors or MacOS.

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
        self.apt_get_update()
        self.apt_get_install_ghostscript()
        self.download_and_install_mafft()

    def apt_get_update(self):
        """Update apt-get package lists"""
        return subprocess.run(["apt-get", "update"])

    def apt_get_install_ghostscript(self):
        """Install ghostscript using apt-get"""
        return subprocess.run(["apt-get", "install", "-y", "ghostscript"])

    def download_and_install_mafft(self, mafft_file="mafft_7.427-1_amd64.deb"):
        """MAFFT exists as a Debian package"""
        deb_filepath = "/tmp/" + mafft_file
        # download the package from mafft.cbrc.jp
        if not os.path.isfile(deb_filepath):
            url = "https://mafft.cbrc.jp/alignment/software/" + mafft_file
            print("Downloading %s..." % url)
            r = requests.get(url, allow_redirects=True)
            open(deb_filepath, "wb").write(r.content)
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
