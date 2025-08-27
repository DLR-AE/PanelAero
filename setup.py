"""
Install Panel Aero with core dependencies via:
- pip install -e <local_repo_path>
"""

from setuptools import setup, find_packages


def my_setup():
    setup(name='PanelAero',
          version='2025.08',
          description='An Implementation of the Vortex Lattice and the Doublet Lattice Method.',
          long_description=open('README.md', encoding='utf8').read(),
          long_description_content_type='text/markdown',
          url='https://github.com/DLR-AE/PanelAero',
          author='Arne Voß',
          author_email='arne.voss@dlr.de',
          license='BSD 3-Clause License',
          packages=find_packages(),
          python_requires='>=3.7',
          install_requires=['numpy'],
          extras_require={'test': ['pytest',
                                   'pytest-cov',
                                   ],
                          'tutorials': ['jupyter',
                                        'jupyter-book',
                                        'matplotlib',
                                        ]
                          },
          )


if __name__ == '__main__':
    my_setup()
