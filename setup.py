from setuptools import setup

setup(
    name='CAVA',
    version='2.0.0',
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
        'bin/ensembl_db',
        'bin/dbSNPDB.py',
        'bin/dbsnp_db'
    ],
    zip_safe=False,
    include_package_data=True
)
