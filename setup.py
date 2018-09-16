try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Process data for BioTree',
    'author': 'Guanglai Li',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'My email.',
    'version': '0.1',
    'install_requires': ['nose', 'pandas', 'numpy', 'matplotlib', 'scipy',
                         'xlrd'],
    'packages': ['biotree'],
    'scripts': [],
    'name': 'biotree',
    'include_package_data': True
}

setup(**config)
