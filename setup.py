from setuptools import setup, find_packages
import os

def read_requirements():
    """Read the requirements.txt file and return a list of dependencies"""
    with open('requirements.txt') as req:
        return req.read().splitlines()

setup(
    name='FastGEMF',
    version='0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=read_requirements(),
    author='Mohammad Hossein Samaei, Faryad Darabi Sahneh, Caterina Scoglio',
    author_email='msamaei@ksu.edu',
    description='Fast Scalabe Spread Process Simulator Over Large Scale Complex Multi-Layer Networks',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/KSUNetSE/fastgemf',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  

)