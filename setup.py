#!/usr/bin/env python
from setuptools import setup, find_packages
import os

scripts = ['bin/exseek']

#with open('requirements.txt') as fin:
#    requirements = [line.strip() for line in fin]
requirements = []

setup(
    name='exseek-pipeline',
    version='1.0.0',
    scripts=scripts,
    packages=['exseek'],
    package_data={'exseek': ['config/*.yaml', 'scripts/*', 'singularity/*', 'templates/*', 'snakefiles/*', 'snakefiles/rules/*', 'snakefiles/scripts/*']},
    include_package_data=True,
    install_requires=requirements,
    url='https://github.com/lulab/exseek',
    description='exSEEK - a pipeline for analysis of exRNA sequencing data',
    zip_safe=False,
    data_files=[('', ['LICENSE', 'requirements.txt'])]
)
