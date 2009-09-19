#!/bin/bash

mkdir run ; rm -rf run/* ; cp grmhd run ; cd run
ln -s ../stellar1.txt .
ln -s ../eosnew.dat . ; ln -s ../eosnew.head . ; ln -s ../eosdegennew.dat .
ln -s ../eosextranew.dat . ; ln -s ../eosextranew.head . ; ln -s ../eosextradegennew.dat .
ln -s ../eossimplenew.dat . ; ln -s ../eossimplenew.head . ; ln -s ../eossimpledegennew.dat .
ln -s ../eossimpleextranew.dat . ; ln -s ../eossimpleextranew.head . ; ln -s ../eossimpleextradegennew.dat .

# for no MPI or OpenMP:
nohup ./grmhd &


