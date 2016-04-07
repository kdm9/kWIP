#########
# SETUP #
#########
set -e
source /hbb_exe_gc_hardened/activate
set -x

export KWIP_VER=$1

cd $(mktemp -d)

########
# DEPS #
########

# Install libbz2
tar xf /io/deps/libbz2.tar*
pushd bzip2*/
make CFLAGS=-fPIC
make install
popd
rm -rf boost*/


#########
# BUILD #
#########

prefix=/io/kwip_${KWIP_VER}_amd64
mkdir -p $prefix

builddir=$(mktemp -d)
pushd $builddir

cmake /io -DCMAKE_INSTALL_PREFIX=$prefix
make -j4
make test
make install
popd

chown 1000:1000 -R $prefix
