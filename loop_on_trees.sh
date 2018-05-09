#!/bin/bash

file="../../../PbPb_2015_TREE/run_list.txt"
while IFS= read line
do
	tree="../../../PbPb_2015_TREE/Tree_$line.root"
	if [ -f "$tree" ]; then
	root -b -q create_tree_from_tree.C++\($line\)
	else
	echo "the file Tree_$line.root does not exist"
	fi
done <"$file"
