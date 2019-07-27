set -xe

# setup
test -d $HOME/miniconda/bin || (wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -u -p $HOME/miniconda)
export PATH="$HOME/miniconda/bin:$PATH"
source $HOME/miniconda/etc/profile.d/conda.sh
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda init

# prepare env
conda create -q -n test-environment numpy cython pytables
conda activate test-environment

# build & run tests
python3 setup.py install
python3 setup.py test

kwipy-count --help
