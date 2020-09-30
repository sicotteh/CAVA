from setuptools import setup
with open('VERSION') as version_file:
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
        'bin/cava',
        'bin/EnsemblDB.py',
        'bin/dbSNPDB.py',
        'bin/dbsnp_db'
    ],
    install_requires = requires,
    zip_safe=False,
    include_package_data=True
)
