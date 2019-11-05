# PyPack manual

## Installing PyPack
PyPack is developed in a linux environment, so if a linux machine is available, this would probably be the best candidate for installing it and running it on. After cloning PyPack to a local directory, `pypack-master/`, you will need to install Plantri. To do this download the gzipped tar file <a href="http://users.cecs.anu.edu.au/~bdm/plantri/plantri50.tar.gz">plantri50.tar.gz</a>, from <a href="http://users.cecs.anu.edu.au/~bdm/plantri/">http://users.cecs.anu.edu.au/~bdm/plantri/</a> into `pypack-master/`. Then open up a bash shell and execute the script `setup.sh` (either by typing `bash setup.sh`, or by first changing the access control bits to allow for execution by typing `chmod 764 setup.sh`, and then running it directly with `./setup.sh`). This will extract all archives and compile Plantri.

## Running PyPack
To enumerate the polyhedra, choose a number N for which all polyhedra having at most N ideal vertices will be enumerated. Then, in a bash shell, navigate to `pypack-master/python/circle_packings/`, and run `enumerate.sh N`. Again, you will probably need to either `chmod` the file permissions, or execute indirectly using the `bash` command. Note that you will need Python and various Python libraries, including Numpy, Matplotlib, Scipy, and cPickle to do this. This will run Plantri N times, and then execute `pack.py` to compute the geometric structures and volumes of the polyhedra from the graphs. This information will be `cPickle`d into the `polyhedra` directory each time another polyhedron is found so that fatal termination of the execution will not lose all of the data.

## Viewing the output
The output polyhedra stored in the `polyhedra` directory and are numbered in the order that they appear in the output from Plantri. To view the pair of circle packings and volume for the n-th polyhedron, navigate to `pypack-master/python/circle_packings/` and run `python read.py n`.
