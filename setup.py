# This setup script uses setuptools to package the geopocket module.
# It specifies the package name, version, and the directory where the source code is located.   
# It also defines an entry point for a console script called 'geopocket' that runs the main function in the predictor module.

from setuptools import setup, find_packages

setup(
    name='geopocket',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'geopocket = geopocket.predictor:main',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'scikit-learn'
    ],
    python_requires='>=3.8',
)
