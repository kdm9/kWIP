
set -xe

cd src
make test
./test
exit 0

test -f miniconda.sh || wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
test -d miniconda || bash miniconda.sh -b -p ./miniconda
export PATH="$PWD/miniconda/bin:$PATH"
conda config --set always_yes true --set changeps1 false
conda update -q conda
conda create -q -n test-environment numpy cython pytables
source activate test-environment

python3 setup.py install
python3 setup.py test

