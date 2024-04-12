# This files describe how GraphAligner handles nodeids from the input graph

This is to make sure that we understand how the seed index links to nodeids and offsets.
Also, to be able to visualize directly in the graph any region where we would observe weird seeds.


## graph loading

Done in `Àligner.cpp` :
`auto alignmentGraph = getGraph(params.graphFile, &memseeder, params);`

This function determines graph format via its extension (.vg or .gfa).
Then, uses either :

```
AlignementGraph result = DirectedGraph::BuildFromVG(graph);

or 

GfaGraph graph = GfaGraph::LoadFromFile(graphFile);
AlignmentGraph result = DirectedGraph::BuildFromGFA(graph);
```

## let's dive in to the GFA case

* Function `GfaGraph::LoadFromFile(filename)`

This function shows a lot of checks, such are refusing empty nodes (sequence="*"), verifying that there is not nodeid duplicates...
Some interesting case that are rejected: "Unspecified edge overlaps (*) are not supported".

  1. loads S and L lines
  2. for S lines, uses `size_t id = getNameId(nameMapping, idstr, result.nodes, result.originalNodeName);`
  3. in getNameId(): `originalNodeName.emplace_back(name);`
  4. ```
     GfaGraph result;
     // S lines
     //In C++ they got the "bright" idea of overloading the rightshift and leftshift operators with streams to represent serialization/deserialization !
     std::stringstream sstr {line};
     std::string idstr;
     std::string dummy;
     std::string seq;
     sstr >> dummy;
     assert(dummy == "S");
     sstr >> idstr;
     size_t id = getNameId(nameMapping, idstr, result.nodes, result.originalNodeName);
     sstr >> seq;
     result.nodes[id] = seq;
     // L lines
     result.edges.emplace_back(frompos, topos, overlap);
     result.edges.emplace_back(topos.Reverse(), frompos.Reverse(), overlap);
     ```
   5. IMPORTANT: id :size_t is the internal id in GraphAligner, which is just a node counter and is mapped to originalNodeName
      For each ne node, `nodeMapping[idstr]=nodeMapping.size()` and `originalNodeName.emplace_back(name)`  <<== the vector position defines the node id
   6. To get original nodeid from GraphAligner internal id use : `GfaGraph::OriginalNodeName(int internalNodeId)`
   7. To get internal id from original nodeid : `GfaGraph graph.originalNodeName[internalNodeId]`
   
* Function `DirectedGraph::BuildFromGFA(graph)`

This one iterates on `GfaGraph graph.edges` and `GfaGraph graph.nodes`. In particuler it calls for `AlignmentGraph result.

  1. internal nodeid from GfaGraph is discarded, original nodeid is passed :
  ```
  for (size_t i = 0; i < graph.nodes.size(); i++)
  {
      std::string name = graph.OriginalNodeName(i);
      result.AddNode((size_t)(i * 2), graph.nodes[i], name, false, breakpointsFw);
  }
  ```
  2. Inside this function :
     `void AlignmentGraph::AddNode(size_t bigraphNodeId, const DNAString& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints)`
     * original node name is copied: `originalNodeName.push_back(name);`
     * and because it is node passed, `GfaGraph` internal nodeid is lost
  3. From `AlignmentGraph`, GFA nodeid should consequently be similar to the result of `AlignmentGraph::BigraphNodeName(size_t bigraphNodeId)`, which just `return originalNodeName[bigraphNodeId];`
  

  
It's a simple temporary wrapper where node / edges are vectors.

```
class GfaGraph
{
public:
	GfaGraph();
	static GfaGraph LoadFromFile(std::string filename);
	static GfaGraph LoadFromStream(std::istream& stream);
	std::string OriginalNodeName(int nodeId) const;
	size_t totalBp() const;
	std::vector<DNAString> nodes;
	std::vector<std::tuple<NodePos, NodePos, size_t>> edges;
	std::vector<std::string> originalNodeName;
private:
};
```


* Function `DirectedGraph::BuildFromGFA(graph)`
