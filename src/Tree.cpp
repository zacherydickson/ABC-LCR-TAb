#include <algorithm>
#include <cstdio>
#include <cstring>
#include "Tree.hpp"

/*############################################################################*/
/*######## TREE NODE #########################################################*/
/*############################################################################*/

/*# CON-/DESTRUCTORS #########################################################*/

//When destroying a node, also destroy all of its children
TreeNode::~TreeNode(){
    for(TreeNode * child : this->children){
        delete child;
    }
    children.clear();
}


/*# METHODS #################################################################*/

//Input string ias assumed to be a cleaned newick, without outer brackets:
//  ie. A,B:6,(C,D:10) instead of (A,B:6,(C,D:10));
void TreeNode::add_children_from_newick(std::string str, unsigned int & nextID){
    while(str.size()){
        //Create the new node and its connections
        TreeNode * child = new TreeNode;
        child->id = id;
        child->parent = this;
        this->children.push_back(child);
        //Extract the string for a single child
        size_t commaPos = std::string::npos;
        int layer = 0;
        int i = 0;
        while(i < str.size() && commaPos == std::string::npos){
            switch(str[i]){
                case '(':
                    layer++;
                    break;
                case ')':
                    layer--;
                    break;
                case ',':
                    if(layer == 0){
                        commaPos = i;
                    }
                    break;
            }
            i++;
        }
        std::string childStr = str.substr(0,commaPos);
        str.erase(0,childStr.size()+1);
        //Child string is in the format (\([A-Za-z,:0-9]\))?[A-Za-z0-9]*(:[0-9]+.?[0-9]*)
        //The branch length is optional
        //The node ID is optional if it is internal
        //First strip off the grandchild string
        size_t closePos = childStr.rfind(')');
        std::string grandchildStr;
        if(closePos != std::string::npos){
            grandchildStr = childStr.substr(0,closePos+1);
            childStr.erase(0,grandchildStr.size());
        }
        //Check if there is a branch length
        size_t branchPos = childStr.rfind(':');
        if(branchPos != std::string::npos){
            std::string branchStr = childStr.substr(branchPos + 1,childStr.size());
            childStr.erase(branchPos,childStr.size());
            child->height = std::stod(branchStr);
        }
        //All that remains is either an empty string or an ID
        if(childStr.size() == 0){ //Need to provide an id
            char buff[100];
            std::sprintf(buff,"NODE_%05d",nextID++);
            child->id = buff;
        } else {
            child->id = childStr;
        }
        //If there is a grandchild string, recurse after stripping leading and trailing
        //parenthesis
        if(grandchildStr.size() > 2){
            child->add_children_from_newick(grandchildStr.substr(1,grandchildStr.size()-2),nextID);
        }
    }
}

cvTreeNode TreeNode::get_all_descendents() const {
    cvTreeNode result;
    for(TreeNode* child : this->children){
        result.push_back(child);
        cvTreeNode cRes = child->get_all_descendents();
        result.insert(result.end(),cRes.begin(),cRes.end());
    }
    return result;
}

cvTreeNode TreeNode::get_descendent_leaves() const {
    cvTreeNode result;
    for(TreeNode* child : this->children){
        if(child->children.size() == 0){
            result.push_back(child);
            continue;
        }
        cvTreeNode cRes = child->get_descendent_leaves();
        result.insert(result.begin(),cRes.begin(),cRes.end());
    }
    return result;
}

std::ostream& TreeNode::output(std::ostream& os) const{
    os << this->id << " @ " << this->height << "\n";
    for(TreeNode* child : this->children){
        child->output(os);
    }
    return os;
}

/*############################################################################*/
/*######## TREE CONTAINER ####################################################*/
/*############################################################################*/

/*# CON-/DESTRUCTORS #########################################################*/

//Given a newick formatted string, parses the string and constructs the
// appropriate Tree
Tree::Tree(std::string newick){
    //Remove everything after a semicolon if present
    newick = newick.substr(0,newick.find(';'));
    //Remove all whitespace characters if present
    char chars [] = " \t\n\r";
    for(unsigned int i = 0; i < strlen(chars); ++i){
        newick.erase(std::remove(newick.begin(), newick.end(), chars[i]), newick.end());
    }
    //Check if there is a label for the root node
    size_t closePos = newick.rfind(')');
    if(closePos < newick.size() - 1){
        this->mRoot.id = newick.substr(closePos+1,newick.size());
        newick.erase(closePos+1,newick.size());
    } else {
        this->mRoot.id = "ROOT";
    }
    if(newick.size() > 2){
        mRoot.add_children_from_newick(newick.substr(1,newick.size()-2),mNextID);
    }
}


/*# METHODS ##################################################################*/

cvTreeNode Tree::get_nodes() const {
    cvTreeNode result = this->mRoot.get_all_descendents();
    result.insert(result.begin(), &this->mRoot);
    return result;
}


std::vector<std::string> Tree::get_leafLabels() const{
    std::vector<std::string> labels;
    for(const TreeNode* node : this->get_leaves()){
        labels.push_back(node->id);
    }
    return labels;
}
