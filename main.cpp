#include <cstdio>
#include <iostream>
#include <iterator>
#include <cstdint>

#include <tomahawk/tomahawk_output_reader.h>
#include "graph.h"


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
	if(graph.build(two_reader) == false){
		std::cerr << "failed" << std::endl;
		return false;
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

	std::cerr << "Linearize" << std::endl;
	if(graph.linearlizeGridMatrix(true) == false){
		std::cerr << "failed linearlization" << std::endl;
		return(1);
	}

	std::cerr << "Cycling over: [" << graph.n_grid_size_ << "," << graph.n_grid_size_ << "] matrix..." << std::endl;
	for(auto it = graph.edges_.begin(); it != graph.edges_.end(); ++it){
		it->computeStatistics();
		javelin::SummaryStatistics& largest_score = it->getLargestScore();

		std::cout <<
				 it->node_in_->parent_contig_->name << "\t" <<
				 it->node_out_->parent_contig_->name << "\t" <<
			(U64)it->left_AA_.n_total << "\t" <<
			(U64)it->left_AB_.n_total << "\t" <<
			(U64)it->left_BA_.n_total << "\t" <<
			(U64)it->left_BB_.n_total << "\t" <<

			(U64)it->right_AA_.n_total << "\t" <<
			(U64)it->right_AB_.n_total << "\t" <<
			(U64)it->right_BA_.n_total << "\t" <<
			(U64)it->right_BB_.n_total << "\t" <<
			     largest_score.mean    << "\t" <<
				 largest_score.total   << "\t" <<
				 largest_score.n_total << "\t" <<
				 largest_score.standard_deviation << std::endl;
	}

	return(0);


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
