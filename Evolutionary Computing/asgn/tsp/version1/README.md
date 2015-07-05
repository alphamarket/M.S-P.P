# Genetic approach for solving TSP problem

Compilation
---
```BASH
cd /path/to/project/
mkdir .build
cd .build
clear; cmake .. && make && ./main [OPTIONS]
```

Accepted Options
---
```
Allowed options:
  --help                                Produce help message
  -d [ --data-file ] arg (=data.dat)    The data input file
  -c [ --config-file ] arg (=config.json)
                                        A json formatted config file
  --cache-file arg                      The cache file to cache loaded data
```
> **Data File** For standard data file see `data/*.tsp`, see [this](http://www.math.uwaterloo.ca/tsp/index.html) for more data

> **Cache File** You can specify the cache file's address, which will be used when after loaded data from `data-file`
to cache loaded data for next uses. Or use the cache file to load previously cached data.
