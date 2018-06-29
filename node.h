#ifndef NODE_H_
#define NODE_H_

namespace javelin{

/**<
 * Node takes edges as references (pointers). Stores edges
 * in a STL-like container (`EdgeReferenceContainer`)`
 */
struct Node{
public:
	typedef Node                    self_type;
	typedef EdgeReferenceContainer  value_type;
    typedef value_type&             reference;
    typedef const value_type&       const_reference;
    typedef value_type*             pointer;
    typedef const value_type*       const_pointer;
    typedef std::ptrdiff_t          difference_type;
    typedef std::size_t             size_type;
    typedef Edge                    edge_type;
    typedef tomahawk::totempole::HeaderContig contig_type;

public:
    Node() : id_(0), parent_contig_(nullptr){}
    Node(const uint32_t id, contig_type* contig) : id_(id), parent_contig_(contig){}
    Node(const self_type& other) :
    	id_(other.id_),
		edges_in_(other.edges_in_),
		edges_out_(other.edges_out_),
		parent_contig_(other.parent_contig_)
    {
    }

    ~Node(){
    	// do not delete parent
    }

	// Capacity
    inline const bool empty(void) const{ return(this->inDegree() == 0 && this->outDegree() == 0); }
    inline const size_type& inDegree(void) const{ return(this->edges_in_.size()); }
    inline const size_type& outDegree(void) const{ return(this->edges_out_.size()); }
    inline const void addEdgeIn(edge_type& edge){ this->edges_in_ += edge; }
    inline const void addEdgeOut(edge_type& edge){ this->edges_out_ += edge; }

    // Accessor
    inline contig_type* getParentContig(void){ return(this->parent_contig_); }
    inline const contig_type* getParentContig(void) const{ return(this->parent_contig_); }

public:
    uint32_t     id_;
	value_type   edges_in_;
	value_type   edges_out_;
	contig_type* parent_contig_;
};

// friend
void Edge::operator+=(const twk_entry_type& twk_entry){
	// Error because Node is defined after this
	const U32& parentA_bp = this->node_in_->parent_contig_->n_bases;
	const U32& parentB_bp = this->node_out_->parent_contig_->n_bases;

	// Assuming configuration AA
	// B |----------* A     A *----------| B
	U32 AA_pos1 = parentA_bp - twk_entry.Aposition;
	U32 AA_pos2 = twk_entry.Bposition;

	// Update only left and right halves given the target contig lengths
	if(AA_pos1 < parentA_bp/2) this->left_AA_ += twk_entry.R2;
	else this->right_AA_ += twk_entry.R2;
	if(AA_pos2 < parentB_bp/2) this->left_AA_ += twk_entry.R2;
	else this->right_AA_ += twk_entry.R2;

	// Assuming configuration AB
	// B |----------* A     B |----------* A
	U32 AB_pos1 = parentA_bp - twk_entry.Aposition;
	U32 AB_pos2 = parentB_bp - twk_entry.Bposition;

	if(AB_pos1 < parentA_bp/2) this->left_AB_ += twk_entry.R2;
	else this->right_AB_ += twk_entry.R2;
	if(AB_pos2 < parentB_bp/2) this->left_AB_ += twk_entry.R2;
	else this->right_AB_ += twk_entry.R2;

	// Assuming configuration BA (as given)
	// A *----------| B     A *----------| B
	U32 BA_pos1 = twk_entry.Aposition;
	U32 BA_pos2 = twk_entry.Bposition;

	if(BA_pos1 < parentA_bp/2) this->left_BA_ += twk_entry.R2;
	else this->right_BA_ += twk_entry.R2;
	if(BA_pos2 < parentB_bp/2) this->left_BA_ += twk_entry.R2;
	else this->right_BA_ += twk_entry.R2;

	// Assuming configuration BB
	// A *----------| B     B |----------* A
	U32 BB_pos1 = twk_entry.Aposition;
	U32 BB_pos2 = parentB_bp - twk_entry.Bposition;

	if(BB_pos1 < parentA_bp/2) this->left_BB_ += twk_entry.R2;
	else this->right_BB_ += twk_entry.R2;
	if(BB_pos2 < parentB_bp/2) this->left_BB_ += twk_entry.R2;
	else this->right_BB_ += twk_entry.R2;

	// Debug
	//std::cerr << parentA_bp << ", " << parentB_bp << "\t" << AA_pos1 << "," << AA_pos2 << "\t" << AB_pos1 << "," << AB_pos2 << "\t" << BA_pos1 << "," << BA_pos2 << "\t" << BB_pos1 << "," << BB_pos2 << std::endl;
}

}



#endif /* NODE_H_ */
