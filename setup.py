import os
from os.path import relpath, join
from setuptools import setup
import versioneer


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


setup(
    name="peleffybenchmarktools",
    author="Martí Municoy",
    author_email="martimunicoy@gmail.com",
    description=("Benchmark tools for peleffy"),
    license="MIT",
    keywords="molecular mechanics, forcefield, potential energy",
    url="https://github.com/martimunicoy/peleffy-benchmarks",
    packages=[
        'peleffybenchmarktools',
        'peleffybenchmarktools/dihedrals'
    ],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Natural Language :: English",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"
    ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
