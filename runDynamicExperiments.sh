#!/bin/bash

exec >dynamic.txt 2>&1
set -x
set -e

make clean
make parallel

echo "Dynamic Scheduling varying thread nums"
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=1
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=2
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=4
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=8
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=16
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=32
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
export OMP_NUM_THREADS=56
numactl --cpunodebind=0 --preferred=0 ./bin/tap 1e-8 1 net/ChicagoRegional_net.txt  net/ChicagoRegional_trips.txt
