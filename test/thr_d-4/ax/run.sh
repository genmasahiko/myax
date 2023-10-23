#!/bin/zsh

for i (disp-*.in) pw.x < $i > $i:r.out
