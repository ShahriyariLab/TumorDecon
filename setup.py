import setuptools

with open("README.md", "r") as fh:
    readme = fh.read()

setuptools.setup(
    name="TumorDecon",
    version="0.2.2",
    author="ShahriyariLab",
    author_email="lshahriyari@umass.edu",
    description="Deconvolution Methods for Digital Cytometry",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/ShahriyariLab/TumorDecon",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    license="MIT License",
    include_package_data=True,
    dependency_links=["git+https://github.com/kristyhoran/singscore"],
    install_requires=[
        "scikit-learn>0.20.3",
        "numpy>1.18.2",
        "pandas>1.0.3",
        "matplotlib>2.2.2",
        "scipy>1.2.1",
        "requests>2.22.0",
        "beautifulsoup4>4.6.0",
        "multiprocess>=0.70.9",
        "seaborn>0.10.0",
        "wget>=3.2",
        "combat>=0.1.7",
        "statsmodels>=0.9.0"
    ]
)
