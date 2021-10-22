from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='cgmbrush',
    version='0.1.0',
    description='CGM Brush, a fast python library for painting baryons in 2D',
    long_description=readme,
    author='Ian Williams',
    author_email='ianw89@live.com',
    requires=['numpy','scipy','pandas','jupyter'],
    url='https://github.com/ianw89/cgm-brush',
    license=license,
    packages=['cgmbrush',
            'cgmbrush.plots',
            ]
)