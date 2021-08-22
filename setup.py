# File: setup.py
# Date: 6-Oct-2018
#
# Update:
#
import re

from setuptools import find_packages
from setuptools import setup
import sys
if sys.version_info[0] < 3:
    from io import open as open

packages = []
thisPackage = "wwpdb.utils.oe_util"

with open("wwpdb/utils/oe_util/__init__.py", "r", encoding="utf-8") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")


setup(
    name=thisPackage,
    version=version,
    description="wwPDB workflow engine utils",
    long_description="See:  README.md",
    author="Ezra Peisach",
    author_email="ezra.peisach@rcsb.org",
    url="https://github.com/rcsb/py-wwpdb_apps_wf_engine_utils",
    #
    license="Apache 2.0",
    classifiers=[
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    entry_points={"console_scripts": []},
    #
    # Also depends on 'openeye.oechem' but cannot install by pypi
    install_requires=[
        "mmcif.utils",
        "wwpdb.utils.cc_dict_util ~= 0.2",
        "OpenEye-toolkits-python2.7-ucs2-linux-x64 == 2016.6.1; python_version == '2.7' and sys_platform=='linux2'",
        "OpenEye-toolkits ; python_version >= '3'",
    ],
    packages=find_packages(exclude=["wwpdb.utils.tests_oe_util"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    # These basic tests require no database services -
    test_suite="wwpdb.utils.tests_oe_util",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={
        "dev": ["check-manifest"],
        "test": ["coverage"],
    },
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
