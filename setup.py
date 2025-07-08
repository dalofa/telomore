from setuptools import setup, find_packages

setup(
    name='Telomore',
    version='0.4',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'telomore': ['bash_scripts/*.sh'],
    },
    install_requires=[
        'Bio==1.8.0',
        'biopython==1.80',
        'GitPython==3.1.44',
        'pysam==0.23.3',
    ],
    py_modules=['telomore'],
    author='David Faurdal',
    author_email='dalofa@biosustain.dtu.dk',
    description='A brief description of your package',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/dalofa/telomore',
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'telomore=telomore.telomore:entrypoint',
        ],
    },
)
