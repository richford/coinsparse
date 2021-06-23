"""A simple command line tool to interact with phenotypic data from the Healthy Brain Network."""
import datetime
import os
import os.path as op
import string

from setuptools import setup
from setuptools_scm import get_version


def local_version(version):
    """Patch in a version that can be uploaded to test PyPI."""
    scm_version = get_version()
    if "dev" in scm_version:
        gh_in_int = []
        for char in version.node:
            if char.isdigit():
                gh_in_int.append(str(char))
            else:
                gh_in_int.append(str(string.ascii_letters.find(char)))
        return "".join(gh_in_int)
    else:
        return ""


# Allow installation without git repository, e.g. inside Docker.
if os.path.exists(".git"):
    opts = dict(
        use_scm_version={
            "root": ".",
            "relative_to": __file__,
            "write_to": op.join("coinsparse", "_version.py"),
            "local_scheme": local_version,
        },
    )
else:
    opts = dict(version="0+d" + datetime.date.today().strftime("%Y%m%d"))


if __name__ == "__main__":
    setup(**opts)
