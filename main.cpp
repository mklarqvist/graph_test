#include <cstdio>
#include <iostream>
#include <iterator>
#include <cstdint>

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
public:
	typedef Graph           self_type;
	typedef Node            node_type;
    typedef Edge            edge_type;
    typedef EdgeContainer   edge_container_type;
    typedef EdgeReferenceContainer   edge_reference_container_type;
    typedef NodeContainer   node_container_type;
    typedef std::size_t     size_type;
    typedef tomahawk::totempole::HeaderContig contig_type;
    typedef tomahawk::io::OutputEntry  twk_entry_type;

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

    //
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

public:
    U32 n_grid_size_;
    edge_container_type edges_;
    node_container_type nodes_;
    edge_type** edge_grid_;
};

}

int main(int argc, char** argv){
	if(argc < 2){
		std::cerr << "need parameters" << std::endl;
		return(1);
	}

	std::cerr << std::string(argv[0]) << std::endl;
	std::cerr << std::string(argv[1]) << std::endl;
	std::string input_file = std::string(argv[1]);

	tomahawk::TomahawkOutputReader two_reader;
	if(two_reader.open(input_file) == false){
		std::cerr <<"Failed to open" << std::endl;
		return 1;
	}

	const U32 n_contigs = two_reader.getHeader().getMagic().getNumberContigs();
	if(n_contigs == 0){
		std::cerr << "No contig data! Corrupted..." << std::endl;
		return(1);
	}

	std::cerr << "Adding: " << n_contigs << " nodes..." << std::endl;
	javelin::Graph graph(n_contigs + 10); // max capacity in ctor
	for(U32 i = 0; i < n_contigs; ++i){
		graph.addNode(&two_reader.getHeader().contigs_[i]);
	}
	if(graph.buildGrid(n_contigs) == false){
		std::cerr << "failed to build grid..." << std::endl;
		return(1);
	}

	U64 b_total_uncompressed = 0;
	U64 b_total_compressed   = 0;

	std::cerr << "Cycling over: " << two_reader.getIndex().getContainer().size() << " chunks..." << std::endl;
	for(U32 i = 0; i < two_reader.getIndex().getContainer().size(); ++i){
		b_total_compressed += two_reader.getIndex().getContainer().at(i).byte_offset_end - two_reader.getIndex().getContainer().at(i).byte_offset;
		b_total_uncompressed += two_reader.getIndex().getContainer().at(i).uncompressed_size;
	}
	std::cerr << "Data: " << b_total_compressed << "b -> " << b_total_uncompressed << "b..." << std::endl;

	while(two_reader.nextBlock()){
		tomahawk::containers::OutputContainerReference o = two_reader.getContainerReference();
		for(U32 i = 0; i < o.size(); ++i){
			//std::cerr << o[i] << std::endl;
			if(o[i].AcontigID == o[i].BcontigID) continue;
			graph += o[i];
		}
	}

	std::cerr << "Cycling over: [" << graph.n_grid_size_ << "," << graph.n_grid_size_ << "] matrix..." << std::endl;
	for(U32 i = 0; i < graph.n_grid_size_; ++i){
		for(U32 j = 0; j < graph.n_grid_size_; ++j){
			if(graph.edge_grid_[i][j].hasConnectivity() == false){
				std::cerr << "filter no connectivity" << std::endl;
				continue;
			}

			graph.edge_grid_[i][j].computeStatistics();
			javelin::SummaryStatistics& largest_score = graph.edge_grid_[i][j].getLargestScore();
			std::cout << i << "\t" << j << "\t" <<
					graph.nodes_[i].parent_contig_->name << "\t" << graph.nodes_[j].parent_contig_->name << "\t" <<
					(U64)graph.edge_grid_[i][j].left_AA_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].left_AB_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].left_BA_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].left_BB_.n_total << "\t" <<

					(U64)graph.edge_grid_[i][j].right_AA_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].right_AB_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].right_BA_.n_total << "\t" <<
					(U64)graph.edge_grid_[i][j].right_BB_.n_total << "\t" <<
					largest_score.mean << "\t" << largest_score.total << "\t" << largest_score.n_total << "\t" << largest_score.standard_deviation << std::endl;
		}
	}


	graph.addNode(javelin::Node(1, &two_reader.getHeader().contigs_[1])).addNode(javelin::Node(2, &two_reader.getHeader().contigs_[2])).addNode(javelin::Node(9,&two_reader.getHeader().contigs_[9]));
	//graph.nodes_ += node1;
	//graph.nodes_ += node2;
	javelin::Edge edge (&graph.nodes_[0], &graph.nodes_[1]);
	javelin::Edge edge2(&graph.nodes_[1], &graph.nodes_[0]);

	graph.addEdge(edge).addEdge(edge2).addEdge(javelin::Edge(&graph.nodes_[0], &graph.nodes_[0])).addEdge(edge).addEdge(javelin::Edge());
	std::cerr << "graph: " << graph.sizeNodes() << "," << graph.sizeEdges() << std::endl;

	std::cerr << "In-degree node 0: " << graph.nodes_[0].edges_in_.size() << std::endl;
	std::cerr << "Adding edge to node 0..." << std::endl;
	graph.nodes_[0].addEdgeIn(edge);
	std::cerr << "In-degree node 0: " << graph.nodes_[0].edges_in_.size() << std::endl;
	std::cerr << "Adding 4 edges to node 0..." << std::endl;
	graph.nodes_[0].addEdgeIn(edge2);
	graph.nodes_[0].addEdgeIn(edge2);
	graph.nodes_[0].addEdgeIn(edge2);
	graph.nodes_[0].addEdgeIn(edge);
	std::cerr << "In-degree node 0: " << graph.nodes_[0].edges_in_.size() << std::endl;
	std::cerr << "Forcing resizing of node 0..." << std::endl;
	graph.nodes_[0].edges_in_.resize();
	std::cerr << "In-degree node 0: " << graph.nodes_[0].edges_in_.size() << std::endl;

	std::cerr << &graph.nodes_[0].edges_in_[0] << '\t' << &edge << "\tEdge 0->in->parent: " << graph.nodes_[0].edges_in_[0].node_in_->parent_contig_->name << std::endl;
	std::cerr << &graph.nodes_[0].edges_in_[1] << '\t' << &edge2 << "\tEdge 1->in->parent: " << graph.nodes_[0].edges_in_[1].node_in_->parent_contig_->name  << std::endl;

	for(auto it = graph.nodes_[0].edges_in_.cbegin(); it != graph.nodes_[0].edges_in_.cend(); ++it){
		std::cerr << &(*it) << '\t' << it->node_in_->id_ << "," << it->node_out_->id_ << std::endl;
	}
	std::cerr << graph.nodes_[0].edges_out_.size() << std::endl;
	graph.nodes_.resize(10000);
	graph.edges_.resize(20000);
	std::cerr << "nodes:" << std::endl;
	for(auto it = graph.nodes_.cbegin(); it != graph.nodes_.cend(); ++it){
		std::cerr << &(*it) << '\t' << it->id_ << "\t" << it->parent_contig_->name << std::endl;
	}

	std::cerr << "edges:" << std::endl;
	for(auto it = graph.edges_.cbegin(); it != graph.edges_.cend(); ++it){
		std::cerr << &(*it) << '\t' << it->node_in_ << '\t' << it->node_out_ << '\t' << it->valid() << std::endl;
	}

	std::cerr << &graph.nodes_[0].edges_in_[0] << std::endl;

	return(0);
}
