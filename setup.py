# +
from __future__ import annotations

from setuptools import find_packages, setup

with open("surfsci/version.py") as f:
    version: dict[str, str] = {}
    exec(f.read(), version)
__version__ = version["__version__"]


setup(
    name="surfsci",
    version=__version__,
    author="Amir Hajibabaei",
    author_email="a.hajibabaei.86@gmail.com",
    description="workflows for surface science",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=["numpy", "scipy", "ase"],
    url="https://github.com/amirhajibabaei/SurfSci",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
