#!/usr/bin/python3

from setuptools import setup

config = {
    'name': 'run-cgprna',
    'description': 'Toolkits for RNA-Seq data analysis using cgpRna',
    'author': 'Yaobo Xu',
    'url': 'https://github.com/cancerit/cgpRna/run_cgprna',
    'download_url': '',
    'author_email': 'cgphelp@sanger.ac.uk',
    'version': '0.1.0',
    'python_requires': '>= 3.5',
    'setup_requires': ['pytest'],
    'install_requires': [],
    'packages': ['run_cgprna'],
    'package_data': {'run_cgprna': ['config/*.json']},
    'entry_points': {
        'console_scripts': ['run-cgprna=run_cgprna.command_line:main'],
    },
}

setup(**config)
