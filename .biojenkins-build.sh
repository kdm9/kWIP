
set -xe

wdir=`pwd`

test -d khmer-src || git clone https://github.com/ged-lab/khmer khmer-src

# update and install khmer
pushd khmer-src
	git fetch --all
	git reset --hard origin/master
	rm -rf $wdir/khmer
	mkdir $wdir/khmer
	pushd lib
		make install PREFIX=$wdir/khmer/
	popd
popd


# Build kmerclust
rm -rf build
mkdir build
pushd build
	cmake $wdir -DKHMER_ROOT=$wdir/khmer/
	make
popd

# check that we can follow the instructions provided for installation
mkdir $wdir/not_home
# Process: extract commands, change home directory to above, run
grep '^    ' $wdir/README.md |sed 's,$HOME,$wdir/not_home,g' | bash -xe
rm -rf $wdir/not_home
