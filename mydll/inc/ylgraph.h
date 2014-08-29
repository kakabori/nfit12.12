

class YYGraphNode{
  void *clientData;
  set<YYGraphEdge> edges;
};



class YYGraphEdge{
  void* clientData;
  YYGraphNode *toNode;
};


class YYGraph {
  set<YYGraphNode> nodes;
  hash_map<YYGraphNodePair, void*> maxtrixTab;


  int getLocalTMatrix(YYGraphNode &node,
		      int spread,
		      cmatrix &mat,
		      set<YYGraphNode*> nodeList);

  int updateLocal(cmatrix &mat, vecotr<YYGraphNode*> nodeList);

  int deleteNode(YYGraphNode* nodep);

  int getEdgeValue(YYGraphNode &node1, YYGraphNode &node2);

  int updateRowBlock();

};

void YYGraph::_getLocalNodeListHelper(YYGraphNode &node,
				      int spread,
				      set<YYGraphNode*> &nodeList)
{
  if( spread >= 0 )
    nodeList.insert( &node );

  if ( spread > 0 )
    {
      set<YYGraphEdge>::iterator it;
      
      it = node.edges.begin();
      
      while( it != node.edges.end() )
	{
	  _getLocalNodeListHelper(node, spread-1, nodeList);

	  it++;
	}
    }
};
 

int YYGraph::getLocalNodeList(YYGraphNode &node,
			      int spread,
			      set<YYGraphNode*> &nodeList)
{
  nodeList.clear();

  _getLocalNodeListHelper(node, spread, nodeList);
  return nodeList.size();
}

int YYGraph::getLocalTMatrix(YYGraphNode &node,
			     int spread,
			     cmatrix &mat,
			     set<YYGraphNode*> &nodeList)
{
  int dim;
  set<YYGraphNode*>::iterator it1;
  set<YYGraphNode*>::iterator it2;
  int index1;
  int index2;

  dim = getLocalNodeList( node, spread, nodeList);
  
  mat.resize(dim, dim);

  it1 = nodeList.begin();
  index1 = 0;
  while( it1 != nodeList.end() )
    {
      it2 = nodeList.begin();
      index2 = 0;

      while( it2 != nodeList.end() )
	{
	  mat.value(index1, index2)
	    = mat.value(index2, index1)
	    = getEdgeVaule( **it1, **it2 );
	  it2++;
	}
      it1++;
    }
};
