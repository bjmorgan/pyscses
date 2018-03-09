import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = 'README.md'

config = {
    'name': 'gwpb',
    'author': 'Georgina L. Wellock',
    'packages': ['project']
}

setup(**config)
