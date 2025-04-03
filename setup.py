# setup.py
from setuptools import setup, find_packages

setup(
    name='geopredictor',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'geopredictor = geopredictor.predictor:main',
            #                      ^ paquete     ^ archivo .py sin extensión
            #                                    ^ función main() dentro de predictor.py
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'scikit-learn'
    ],
    python_requires='>=3.8',
)

