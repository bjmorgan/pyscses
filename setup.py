from setuptools import setup
from pyscses import __version__ as VERSION
import re

def remove_img_tags(data):
    p = re.compile(r'<img.*?/>')
    return p.sub('', data)

readme = 'README.md'
long_description = open( readme ).read()
# Removing images seems easier than trying to get them to work on PyPI.
long_description = remove_img_tags( long_description )

config = {
    'description': 'PYthon Space-Charge Site Explicit Solver',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'name': 'pyscses',
    'author': 'Georgina L. Wellock',
    'author_email': 'g.l.wellock@bath.ac.uk',
    'packages': ['pyscses'],
    'url': 'https://github.com/bjmorgan/pyscses',
    'download_url': 'https://github.com/bjmorgan/pyscses/archive/%s.tar.gz' % (VERSION),
    'version': VERSION,
    'install_requires': open( 'requirements.txt' ).read(),
    'license': 'MIT'
}

setup(**config)
