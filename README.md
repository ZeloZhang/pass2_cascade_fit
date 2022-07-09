
1) for c++ compilation run 'make'

2) if analysis is needed from python, create shared library in ./src:
g++ -shared -Wl,-soname,analysis.so -o analysis.so *.o ./models/*.o ./systematics/*.o -lpython2.7 -lboost_python `root-config --libs`
