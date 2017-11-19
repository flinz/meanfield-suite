from setuptools import find_packages, setup

with open('README.md') as f:
    readme = f.read()

install_requires = [
    'numpy==1.13.3',
    'scipy==1.0.0',
    'matplotlib==2.1.0',
    'Brian2==2.1.2',
    'pandas==0.21.0',
]

tests_require = [
    'pytest'
]

setup_requires = [
    'pytest-runner',
    'autopep8'
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
