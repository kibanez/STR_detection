from setuptools import setup, find_packages

setup(
    name='STR_detection',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/genomicsengland/STR_detection',
    install_requires=[
        'matplotlib==2.2.2',
        'numpy==1.14.3',
        'subprocess32==3.5.1',
        'pyparsing==2.2.0',
        'setuptools-scm==2.1.0',
        'backports.functools-lru-cache==1.5',
        'cycler==0.10.0',
        'kiwisolver==1.0.1'
        ],
    license='',
    classifiers=[
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Bioinformatics Data science',
        'Topic :: Software Development :: STRs Bioinformatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7.10',
        'Programming Language :: Python :: 3.6.5',
    ],
    author='kristina ibanez garikano',
    author_email='kristina.ibanez-garikano@genomicsengland.co.uk',
    description='Detection of short tandem repeats from NGS'
)