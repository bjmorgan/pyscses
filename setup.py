import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = 'README.md'

config = {
    'name': 'pyscses',
    'author': 'Georgina L. Wellock',
    'packages': ['pyscses'],
    'install_requires': [ 'scipy', 'sympy', 'numpy' ]
}

setup(**config)
