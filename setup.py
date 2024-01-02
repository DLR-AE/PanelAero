"""
Setup file, currently supports:

- installation via "pip install --user <repo_url>"
- installation via "python setup.py install --user"
"""

from setuptools import setup, find_packages


def my_setup():
    setup(name='Panel-Aero',
          version='2023.04',
          description='An Implementation of the Vortex Lattice and the Doublet Lattice Method.',
          url='https://github.com/DLR-AE/PanelAero',
          author='Arne Voß',
          author_email='arne.voss@dlr.de',
          license='BSD 3-Clause License',
          packages=find_packages(),
          python_requires='>=3.7',
          install_requires=[
              'numpy',
              'pytest',
              'pytest-cov',],
          )


if __name__ == '__main__':
    my_setup()
