import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='thermocouples_reference',
    version='0.20',
    description='Thermocouple emf reference functions',
    long_description=read('README.rst'),
    author='User:Nanite @ wikipedia',
    license='public domain',
    url='https://pypi.python.org/pypi/thermocouples_reference',
    packages=['thermocouples_reference'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Manufacturing',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
    ],
    keywords=[
        'thermocouple',
        'thermometer',
        'temperature',
        'emf',
        'electromotive',
        'thermoelectric',
        'Seebeck',
        'lookup',
        'table',
        'NIST',
        'ASTM',
    ],
    install_requires=[
        'numpy',
    ],
    extras_require = {
        'inverse_lookup':  ['scipy'],
    },
    zip_safe=True)

