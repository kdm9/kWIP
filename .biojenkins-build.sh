
wdir=`pwd`

test -d khmer-src || git clone https://github.com/ged-lab/khmer khmer-src

# update and install khmer
pushd khmer-src
	git fetch --all
	git reset --hard origin/master
	rm -rf $wdir/khmer-target
	mkdir $wdir/khmer-target
	pushd lib
		make install PREFIX=$wdir/khmer-target/
	popd
popd


rm -rf build
mkdir build && cd build

cmake $wdir -DKHMER_ROOT=$wdir/khmer-target/
make

# check that we can follow the instructions provided for installation
mkdir $wdir/not_home
# Process: extract commands, change home directory to above, run
grep '^    ' README.md |sed 's,$HOME,$wdir/not_home,g' | bash -xe
rm -rf $wdir/not_home
