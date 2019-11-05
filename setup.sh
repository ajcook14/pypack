#! /bin/bash

tar -xzvf plantri50.tar.gz
mv plantri50 plantri
rm plantri50.tar.gz
cd plantri/
make plantri
mkdir ../python/circle_packings/polyhedra/
