#!/bin/bash

cwd=$(pwd)

cd $cwd/../

pwd

for folder in */; do

	model=${folder%?} # remove the last character (/)

	for file in $( ls $cwd/*$model* ); do
	
		mv $file $cwd/../$model/

	done

done
