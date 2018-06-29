#ifndef GRAPH_H_
#define GRAPH_H_

#include <tomahawk/tomahawk_output_reader.h>
#include "edge.h"
#include "edge_reference_container.h"
#include "edge_container.h"
#include "node.h"
#include "node_container.h"

namespace javelin{

/**<
 * Graph: edge-driven -> G = (V,E)
 */
class Graph{
private:
	typedef Graph                  self_type;
	typedef Node                   node_type;
	typedef Edge                   edge_type;
	typedef EdgeContainer          edge_container_type;
	typedef EdgeReferenceContainer edge_reference_container_type;
	typedef NodeContainer          node_container_type;
	typedef std::size_t            size_type;
	typedef tomahawk::totempole::HeaderContig contig_type;
	typedef tomahawk::TomahawkOutputReader    two_reader_type;
	typedef tomahawk::io::OutputEntry         twk_entry_type;

public:
    Graph(void) :
    	n_grid_size_(0),
		edge_grid_(nullptr)
	{

	}

    Graph(size_type capacity) :
    	n_grid_size_(0),
    	edges_(capacity),
		nodes_(capacity),
		edge_grid_(nullptr)
    {

    }

    ~Graph(){
    	this->clearEdgeGrid();
    }

    /**<
     * Constructs the [N,N] grid-matrix of edges and N vector of vertices
     * @param two_reader Input Tomahawk TWO reader instance
     * @return           Returns TRUE upon success or FALSE otherwise
     */
    bool build(const two_reader_type& two_reader){
    	const U32 n_contigs = two_reader.getHeader().getMagic().getNumberContigs();
		if(n_contigs == 0){
			std::cerr << "No contig data! Corrupted..." << std::endl;
			return(false);
		}

		std::cerr << "Adding: " << n_contigs << " nodes..." << std::endl;

		for(U32 i = 0; i < n_contigs; ++i){
			this->addNode(&two_reader.getHeader().contigs_[i]);
		}

		// Cleanup potential prior grid
		this->clearEdgeGrid();

		if(this->buildGrid(n_contigs) == false){
			std::cerr << "failed to build grid..." << std::endl;
			return(false);
		}

		return(true);
    }

    /**<
     * Construct a [N,N]-matrix of edges
     * @param n_cells Number of columns/rows
     * @return        Returns TRUE upon success or FALSE otherwise
     */
    bool buildGrid(const U16 n_cells){
    	this->n_grid_size_ = n_cells;
    	if(this->nodes_.size() != n_cells){
    		std::cerr << "mismatch between cell size and number of nodes: " << n_cells << "/" << this->nodes_.size() << std::endl;
    		return false;
    	}

    	this->edge_grid_ = static_cast<edge_type**>(::operator new[](n_cells*sizeof(edge_type)));
    	for(U32 i = 0; i < n_cells; ++i){
    		this->edge_grid_[i] = static_cast<edge_type*>(::operator new[](n_cells*sizeof(edge_type)));
    		for(U32 j = 0; j < n_cells; ++j){
    			new( &this->edge_grid_[i][j] ) edge_type( &this->nodes_[i], &this->nodes_[j] );
    		}
    	}

    	return true;
    }

    /**<
     * Collapsed a [N,N]-matrix of edges into a linear container of edges.
     * @param filter_empty Should empty values be filtered out? Defaults to TRUE
     * @return             Returns TRUE upon success or FALSE othewise
     */
    bool linearlizeGridMatrix(const bool filter_empty = true){
    	for(U32 i = 0; i < this->n_grid_size_; ++i){
			for(U32 j = 0; j < this->n_grid_size_; ++j){
				if(this->edge_grid_[i][j].hasConnectivity() == false && filter_empty)
					continue;

				this->edges_ += this->edge_grid_[i][j];
			}
    	}
    	std::cerr << "edges now: " << this->edges_.size() << " down from " << this->n_grid_size_*this->n_grid_size_ << std::endl;
    	return(true);
    }

    inline void operator+=(const twk_entry_type& twk_entry){ this->edge_grid_[twk_entry.AcontigID][twk_entry.BcontigID] += twk_entry; }

    // Capacity
    inline const size_type sizeNodes(void) const{ return(this->nodes_.size()); }
    inline const size_type sizeEdges(void) const{ return(this->edges_.size()); }
    inline const bool empty(void) const{ return(this->sizeNodes() == 0); }

    // Incremental construction
    inline self_type& addNode(const node_type& node){ this->nodes_ += node; return(*this); }
    inline self_type& addNode(const node_type* const node){ this->nodes_ += *node; return(*this); }
    inline self_type& addNode(contig_type* contig){ this->nodes_ += node_type(this->nodes_.size(), contig); return(*this); }
    inline self_type& addEdge(const edge_type& edge){ this->edges_ += edge; return(*this); }
    inline self_type& addEdge(const edge_type* const edge){ this->edges_ += *edge; return(*this); }

    /**<
     * Destroys the current [N,N]-edge grid matrix if available
     */
    void clearEdgeGrid(void){
    	if(edge_grid_ != nullptr){
    		for(U32 i = 0; i < this->n_grid_size_ ; ++i){
    			for(U32 j = 0; j < this->n_grid_size_ ; ++j){
					((&this->edge_grid_[i][j])->~Edge());
				}
				::operator delete[](static_cast<void*>(this->edge_grid_[i]));
			}
			::operator delete[](static_cast<void*>(this->edge_grid_));
		}
    }

public:
    U32 n_grid_size_;
    edge_container_type edges_;
    node_container_type nodes_;
    edge_type** edge_grid_;
};

}



#endif /* GRAPH_H_ */
