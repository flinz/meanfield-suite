
from setuptools import find_packages, setup

readme = open('README.md').read()

install_requires = [

]

tests_require = [

]

setup(
    name='meanfield_suite',
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
