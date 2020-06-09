import setuptools

with open("TumorDecon/README.md", "r") as fh:
    readme = fh.read()

setuptools.setup(
    name="TumorDecon-raronow",
    version="0.0.1",
    author="ShahriyariLab",
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
    license='LICENSE',
    include_package_data=True
)
