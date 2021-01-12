import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="varif",
    version="0.0.3",
    author="Marc-Antoine Guery",
    author_email="marcantoine.guery@gmail.com",
    description="A variant filtering and annotating package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marcguery/varif",
    packages=setuptools.find_packages(),
    scripts=['bin/varif'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)