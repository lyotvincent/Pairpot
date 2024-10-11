#python setup.py sdist bdist_wheel

from setuptools import setup, find_packages
print("packages:",find_packages())

setup(
    name="pairpot",
    version="0.0.3",
    author="Zhihan Ruan, Zhenjie Zhang, and Jian Liu.",
    author_email="rrrzhan@mail.nankai.edu.cn",
    description="Pairpot: a database with real-time lasso-based analysis tailored for paired single-cell and spatial transcriptomics.",
    long_description=open('readme.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/lyotvincent/Pairpot",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        '':['pairpot/resources/Cell_marker_Seq.xlsx'],
        '':['pairpot/resources/PanglaoDB_markers_27_Mar_2020.tsv'],
    },
    install_requires=[
        'h5py',
        'numpy',
        'pandas',
        'anndata',
        'scanpy',
        'pkg_resources',
        'scipy',
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
