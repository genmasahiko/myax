#!/bin/zsh

for i in disp_*; do
	cd $i
	pw.x < in > out
	cd ../
done
