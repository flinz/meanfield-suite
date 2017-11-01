
from setuptools import find_packages, setup

readme = open('README.md').read()

install_requires = [
    'numpy',
    'scipy',
    'matplotlib',
    'brian2'
]

tests_require = [
    'pytest'
]

setup(
    name='meanfield',
    version='0.0.1',
    description='Meanfield suite',
    long_description=readme,
    #author='',
    #author_email='',
    url='https://github.com/flinz/meanfield-suite',
    #license='',
    packages=find_packages(),
    install_requires=install_requires,
    tests_require=tests_require,
)
