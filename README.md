## Introduction
This is the source code of the ICDE 2024 paper "*Scalable Distance Labeling Maintenance and Construction for Dynamic Small-World Networks*" (submitted). Please refer to the paper for the algorithm details.

## Algorithms

The implementation code includes the index construction, query processing, and index update of our *DCT* algorithm. The runnable components include:

1. Tree index strategy: *local tree index*, *global tree index*
1. Core index construction methods: *BPCL*, *PCL*, *PLL*, *PSL*, *GLL*
1. Core index maintenance methods: *PDPLL*, *SDPLL*





## Data
An example graph *GOOGLE* is provided in the directory *data* for your reference. You can run *DCT* on the example graph by using the source path `./data`. Note that you should `unzip -q GOOGLE.zip` to get the processed graph.


## Dependency

1. `g++` and `boost`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
