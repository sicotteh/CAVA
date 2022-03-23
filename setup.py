from setuptools import setup
import os
with open(os.path.join(os.path.dirname(__file__), 'cava', 'VERSION')) as version_file:
    version = version_file.read().strip()

with open('requirements.txt') as requires_file:
    requires = requires_file.read().split('\n')

setup(
    name='CAVA',
    python_requires=">=3.7",
    version=version,
    description='CAVA (Clinical Annotation of VAriants)',
    url='https://github.com/Steven-N-Hart/CAVA',
    author='Steven-N-Hart, Hugues-Sicotte',
    author_email='hart.steven@mayo.edu, sicotte.hugues@mayo.edu',
    license='MIT',
    packages=['cava', 'cava.utils', 'cava.ensembldb'],
    install_requires=requires,
    zip_safe=False,
    include_package_data=True,
    test_suite='tests',
    tests_require=['wget'],
)
