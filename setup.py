"""Python setup.py for depmap_omics_upload package"""
import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("depmap_omics_upload", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
            os.path.join(os.path.dirname(__file__), *paths),
            encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="depmap_omics_upload",
    version=read("depmap_omics_upload", "VERSION"),
    description="depmap_omics_upload created by @5im1z at broadinstitute. Includes loading and uploading functions critical to the DepMap omics pipeline.",
    url="https://github.com/broadinstitute/depmap_omics_upload/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="broadinstitute",
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={
        "console_scripts": ["depmap_omics_upload = depmap_omics_upload.__main__:main"]
    },
    package_data={"depmap_omics_upload": ["depmap_omics_upload/config.json"]},
    include_package_data=True,
    extras_require={"test": read_requirements("requirements-test.txt")},
)
