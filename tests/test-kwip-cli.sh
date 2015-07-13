#!/bin/bash

set -e

cli=bin/kwip

function filesize() {
	wc -c "$1" | awk '{print $1}'
}

tmpdir=${TMPDIR:-/tmp}/$$-$RANDOM-kwip-tests
mkdir $tmpdir
trap "echo rm -rf $tmpdir" EXIT

for unweighted in '' '-U'
do
	tst=empty-default-$kernel-$RANDOM
	out=$tmpdir/${tst}.out
	err=$tmpdir/${tst}.err
	k=$tmpdir/${tst}.kernel
	d=$tmpdir/${tst}.dist
	set -x
	$cli 					\
		-k $k 				\
		-d $d 				\
		-t 1				\
		$unweighted			\
		data/defined-[123].ct 		\
		2>&1 >$out | tee $err
	test $(filesize $k) -gt 0
	test $(filesize $d) -gt 0
	set +x
done
