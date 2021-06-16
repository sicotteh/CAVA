from setuptools import setup
import os
with open(os.path.join(os.path.dirname(__file__), 'bin', 'VERSION')) as version_file:
    version = version_file.read().strip()

with open('requirements.txt') as requires_file:
    requires = requires_file.read().split('\n')

setup(
    name='CAVA',
    version=version,
    description='CAVA (Clinical Annotation of VAriants)',
    url='https://github.com/Steven-N-Hart/CAVA',
    author='Steven-N-Hart',
    author_email='hart.steven@mayo.edu',
    license='MIT',
    packages=['cava_', 'ensembldb'],
    scripts=[
        'bin/CAVA.py',
        'bin/EnsemblDB.py',
        'bin/dbSNPDB.py',
    ],
    install_requires=requires,
    zip_safe=False,
    include_package_data=True
)
