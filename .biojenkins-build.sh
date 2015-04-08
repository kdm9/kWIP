
test -d khmer-src || git clone https://github.com/ged-lab/khmer khmer-src

# update and install khmer
pushd khmer-src
	git fetch --all
	git reset --hard origin/master
	rm -rf ../khmer-target
	mkdir ../khmer-target
	pushd lib
		make install PREFIX=../../khmer-target/
	popd
popd


mkdir build && cd build

cmake .. -DKHMER_ROOT=../khmer-target/
make
