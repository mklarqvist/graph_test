```bash
git clone https://github.com/mklarqvist/graph_test
cd graph_test
git clone --recursive https://github.com/mklarqvist/tomahawk
cd tomahawk
make
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
export LD_LIBRARY_PATH
cd ../Debug
make
```

Usage: all data written to standard out
```
graph <file.two>
```
