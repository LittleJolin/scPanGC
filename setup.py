from setuptools import setup, find_packages

setup(
    name="scPanGC",
    version="0.1.0",
    author="Jolin",
    author_email="I will change it after blind peer review",
    description="A computational framework for extracting Gene Clusters from large-scale single-cell autoimmune atlases.",
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/LittleJolin/scPanGC", 
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
    install_requires=[
        "scanpy>=1.9.0",
        "anndata>=0.8.0",
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "metacells>=0.9.0",
        "scipy>=1.7.0"
    ],
)