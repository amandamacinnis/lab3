#! /bin/bash -u

for file in $(ls -1 data/atlas*.fits)
do
	solve-field --ra 57.291 --dec 24.053 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/hd23778*.fits)
do
	solve-field --ra 57.145 --dec 24.181 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/hd23479*.fits)
do
	solve-field --ra 56.567 --dec 24.19 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/alcyone*.fits)
do
	solve-field --ra 56.871 --dec 24.105 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/celaeno*.fits)
do
	solve-field --ra 56.201 --dec 24.289 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/merope*.fits)
do
	solve-field --ra 56.582 --dec 23.948 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/hd23511*.fits)
do
	solve-field --ra 56.664 --dec 24.103 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

for file in $(ls -1 data/electra*.fits)
do
	solve-field --ra 56.219 --dec 24.113 --radius 0.5 --scale-units degwidth --scale-low 0.25 --scale-high 1.0  ${file}
done

