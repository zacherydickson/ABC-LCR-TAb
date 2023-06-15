#ifndef _TREE_HPP
#define _TREE_HPP

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct TreeNode;

typedef std::vector<const TreeNode*> cvTreeNode; 
typedef std::vector<TreeNode*> vTreeNode; 

struct TreeNode {
    ~TreeNode();
    double height = 0;
    std::string id = "";
    TreeNode* parent = NULL;
    vTreeNode children;
    void add_children_from_newick(std::string str, unsigned int & nextID);
    cvTreeNode get_all_descendents() const;
    cvTreeNode get_descendent_leaves() const;
    std::ostream& output(std::ostream& os) const;
};


class Tree {
    //Con-/Destruction
    Tree() = delete;
    public:
        Tree(std::string newick);
//        Tree(const TreeNode & root);
//        Tree(const Tree & other) : Tree(other.mRoot) {}
        ~Tree() {}
    //Members
    private:
        unsigned int mNextID = 0;
        TreeNode mRoot;
//        std::map<std::string,TreeNode*> mNodeMap;
    //Accessors
    public:
        const TreeNode * get_root() const {return &mRoot;}
//        const TreeNode * get_node(const std::string & id) const;
        cvTreeNode get_nodes() const;
        cvTreeNode get_leaves() const {return this->mRoot.get_descendent_leaves();}
        std::vector<std::string> get_leafLabels() const;
        size_t getNLeaves() const {return this->mRoot.get_descendent_leaves().size();}
    //Methods
    public:
        std::ostream& output(std::ostream& os) const {return this->mRoot.output(os);}
    //Operators
};

inline std::ostream& operator<<(std::ostream& os, const Tree & obj){return obj.output(os);}

template<typename T>
struct SBasicNode {
    double height = 0;
    int parent = -1;
    std::string label;
    T value;
};

template<typename T>
struct SDepthFirstAccessNodeVector {
    //Cons-/Destruction
        SDepthFirstAccessNodeVector(const Tree & tree){
            cvTreeNode vTNodes = tree.get_nodes();
            size_t nNodes = vTNodes.size();
            //Construct a vector of ParentIndexes
            //  For node i, insert nChildren copies of i into the parent vector
            //  Ex. 0 -> 1,3; 1->2; 3->4,5,6
            //  Node Vector: 0,1,2,3,4,5,6
            //  -1 -> -1,0,0 -> -1,0,1,0,3,3,3
            this->vNodes.push_back(SBasicNode<T>());
            for(int i = 0; i < nNodes; i++){
                const TreeNode* & node = vTNodes[i];
                vNodes[i].label = node->id;
                vNodes[i].height = node->height;
                for(const TreeNode* child : node->children){
                    SBasicNode<T> iNode;
                    iNode.parent = i;
                    this->vNodes.push_back(iNode);
                }
                idMap[node->id] = i;
            }
        }
    //Members
        std::vector<SBasicNode<T>> vNodes;
        std::map<std::string,int> idMap;
    //Methods
        size_t size(){return this->vNodes.size();}
        SBasicNode<T> & atLabel(std::string label){return this->vNodes[idMap.at(label)];}
        SBasicNode<T> & atIndex(int idx){return this->vNodes[idx];}
};

#endif



