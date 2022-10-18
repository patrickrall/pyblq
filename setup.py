from setuptools import setup

VERSION = "0.0.1"

with open("requirements.txt") as f:
    REQUIREMENTS = f.read().splitlines()

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pybloq",
    version=VERSION,
    description="Bloq Encoding for Qubits",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/patrickrall/pyblq",
    author="Patrick Rall",
    license="MIT License",
    python_requires=">=3.7",
    include_package_data=True,
    install_requires=(REQUIREMENTS)
    )
