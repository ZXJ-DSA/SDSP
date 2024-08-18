## Introduction
This is the implementation code of the ICDE 2024 paper "*Scalable Distance Labeling Maintenance and Construction for Dynamic Small-World Networks*" . Please refer to the paper for the algorithm details.

## Algorithms

The implementation code includes the index construction, query processing, and index update of our *DCT* algorithm. The runnable components include:

1. Tree index strategy: *local tree index*, *global tree index*
1. Core index construction methods: *BPCL*, *PCL*, *PLL*, *PSL*, *GLL*
1. Core index maintenance methods: *PDPLL*, *SDPLL*





## Data
The datasets of this paper are sourced from [http://snap.stanford.edu/data](http://snap.stanford.edu/data), [http://konect.cc/](http://konect.cc/), and [https://law.di.unimi.it/datasets.php](https://law.di.unimi.it/datasets.php). Please refer to the paper for details.

For your reference, an example graph *GOOGLE* is provided in the directory *data*. You can run *DCT* on the example graph by using the source path `./data`. Note that you should `unzip -q GOOGLE.zip` to get the processed graph.


## Dependency

1. `g++` and `boost`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.




## Reference

If you found any of our implementation useful, please cite the appropriate references, listed below:

```
@inproceedings{zhou2024scalable,
  title={Scalable distance labeling maintenance and construction for dynamic small-world networks},
  author={Zhou, Xinjie and Zhang, Mengxuan and Li, Lei and Zhou, Xiaofang},
  booktitle={2024 IEEE 40th International Conference on Data Engineering (ICDE)},
  pages={4573--4585},
  year={2024},
  organization={IEEE}
}
```
