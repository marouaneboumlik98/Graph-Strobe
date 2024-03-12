# Comments on GraphAligner code

This branch contains information related to the indexation and associated data structures in GraphAligner.
It aims to evaluate how feasible mer scheme replacement can be done in GraphAligner (GA) pipeline.

## Sources overview

### dependencies

GA depends on 4 GIT submodules: 
- BBHash: BBHash is a simple library for building minimal perfect hash function. It is designed to handle large scale datasets. The function is just a little bit larger than other state-of-the-art libraries, it takes approximately 3 bits / elements (compared to 2.62 bits/elem for the emphf lib), but construction is faster and does not require additional memory.
- MEMfinder: from the same author, no documentation. Searches for MEMs (minimal exact matches). I can see it is based on a FM-index (BWT & wavelet tree).
- concurrentqueue: from the author: "An industrial-strength lock-free queue for C++, knock-your-socks-off blazing fast performance." See its repo for more.
- parallel-hashmap: A repository aiming to provide a set of excellent hash map implementations, as well as a btree alternative to std::map and std::set. 
- zstr: Enables the use of C++ standard iostreams to access ZLib-compressed streams.

To pull them: `git submodule pull` after cloning this repository

### Important source files:

- `src/AlignerMain.cpp`: 
    - contains the `main()`
    - setup execution parameters which are registered in variable `AlignerParams params`.
    - then call `alignReads(params)` from `src/Aligner.h`
  - `src/Aligner.cpp`:
    - `void alignReads()` : does the top operation: seeds creation or loading via index, then multithreaded alignment + output writing.

## Seeds computations

  1. `alignReads()` in `src/Aligner.cpp` does the job. 
    2. Next block is the seed creation step:
      ```C++
        MEMSeeder* memseeder = nullptr;
        auto alignmentGraph = getGraph(params.graphFile, &memseeder, params);
        ## ignore for now
        DiploidHeuristicSplitter diploidHeuristic;
        if (params.useDiploidHeuristic) { [...] }
        bool loadMinimizerSeeder = params.minimizerSeedDensity != 0;
        MinimizerSeeder* minimizerseeder = nullptr;
        if (loadMinimizerSeeder)
        {
              std::cout << "Build minimizer seeder from the graph" << std::endl;
              minimizerseeder = new MinimizerSeeder(alignmentGraph, params.minimizerLength, params.minimizerWindowSize, params.numThreads, 1.0 - params.minimizerDiscardMostNumerousFraction);
        }
      ```
      - MEMs are created from the graph each time even is minimizers are called for ?
      - NO! It depends on parameters and the format of the graph to load, everything is tested in getGraph().
    3. `getGraph()` will look for `.vg` or `.gfa` in the input graph filename.
      - VG & MEMs: instanciate the MEMSeeder and calls for `DirectedGraph::BuildFromVG()`
      - VG & !MEMs: calls for `DirectedGraph::StreamVGGraphFromFile()`
      - GFA: calls for `DirectedGraph::BuildFromGFA(graph)` (and MEMSeeder if requested)  
         - The structure to model the graph is `class DirectedGraph`, from `src/BigraphToDigraph.h`
         - It does not model VG/GFA paths !
         - If path are desired, just adding a definition there is enough
           But VG loading functions will need to be extended.
           It requires to load batches of node/edges which need to be merged to build the full graph.
           This is similar to what I tested with the python API, see my code there if necessary.
    4. This loads a Bigraph using the temporary struct defined in `src/BigraphToDigraph.h`.
       Next is a conversion to get a Digraph, functions for conversion are in the same file.
       The results is an object `AlignmentGraph`, the graph structure used in the algorithms.
    5. Back to seed creation step. We have now a `AlignmentGraph alignmentGraph` (and maybe MEMSeeder intialized).
       Next is the initialization of the `MinimizerSeeder` object, defined in `src/MinimizerSeeder.h`.
       - happens only if no seed file loaded.
       - skipping the part related to multithreading
       - for 1 thread, the main loop is: 
         ```C++
           while (true)
               {
                 std::string sequence = graph.BigraphNodeSeq(nodeId);
                 iterateMinimizers(sequence, minimizerLength, windowSize, [this, &nodeMinimizerStart, &positionDistributor, &kmerPerBucket, &positionPerBucket, &vecPos, positionSize, thread, nodeId](size_t pos, size_t kmer)
                 {
                   #computes minimizers for each window
                 });
               }
         ```
       - we expected it, as paths are not loaded. Minimizers are computed per node sequence.
       - `iterateMinimizers()` is a callback either to a) `iterateMinimizersSimple()` or b) `iterateMinimizersReal()`. 
         - a) k-mer are computed using `charToInt()` definitions and one bit left shift.
         - a) hash is done via `hash()`, defined in `MinimizerSeeder.cpp` 
         - a) NOTE: there is a IFDEF for osx hash() in `CommonUtils.h`  
         - b) k-mer are computed via bit shift, mask and bitwise OR
         - b) same hash function

## Seeds from a loaded file

These variables appear to contain seeds loaded from a file (not computed):
```C++
const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
std::unordered_map<std::string, std::vector<SeedHit>> seedHits;
```

6. Back to seed creation step.
  - if a seed file is given via parameters, it loads using a lambda iterating on `vg::Alignment`, it feels the map `seedHits`
  - this appears to be also feasible via a GAF file via `loadGafSeeds()`  
  - then it instanciates a `Seeder seeder` object to which is potentially attach our potentially just computed memseeder or minimizerseeder.
  - depending on `seeder.mode` it will use MUM, MEM or minimizers.
7. Then comes the alignments and output creation


## Alignments

- Looking at ``
- It is multithreaded
- There are output conumer depending on output format, for as many classes.
  - fastqThread
  - GAMwriterThread 
  - GAFwriterThread 
  - JSONwriterThread 
  - correctedWriterThread 
  - correctedClippedWriterThread
- the main running Component for mapping is `runComponentMappings()`
  - there lies the while loop that ierates on sequences
  - if seeder.mode!=None, query seeds are generated via this line:
    `std::vector<SeedHit> seeds = seeder.getSeeds(fastq->seq_id, fastq->sequence);`
  - this returns a `std::vector<SeedHit> seeds`
  - `SeedHit` defined in `GraphAlignerWrapper.h`, this is a general container to get:
    ```C++
    int nodeID;
    size_t nodeOffset;
    size_t seqPos;
    size_t matchLen;
    bool reverse;
    size_t rawSeedGoodness;  #<== not sure what is that
    ```
  - `Seeder.getSeeds()` calls to one of the registered seeder, for instance:
    ```C++
    case Mode::Minimizer:
      assert(minimizerSeeder != nullptr);
      return minimizerSeeder->getSeeds(seq, minimizerSeedDensity);
    ```
  - In this function, after a callback to `iterateKmers()`
  - It seems we call a similar mer coding and hashing
    