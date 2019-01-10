from setuptools import setup
from pyscses import __version__ as VERSION

readme = 'README.md'
long_description = open( readme ).read()

config = {
    'description': 'PYthon Space-Charge Site Explicit Solver',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'name': 'pyscses',
    'author': 'Georgina L. Wellock',
    'packages': ['pyscses'],
    'url': 'https://github.com/bjmorgan/pyscses',
    'download_url': 'https://github.com/bjmorgan/pyscses/archive/%s.tar.gz' % (VERSION),
    'version': VERSION,
    'install_requires': open( 'requirements.txt' ).read(),
    'license': 'MIT'
}

setup(**config)
