#ifndef EDGE_H_
#define EDGE_H_

namespace javelin{

//forward declare
struct Node;

// Connectivity
// AA, AB, BA, BB
// r2 bins -> counts, mean, min, max, sd
// inferred length distribution
struct SummaryStatistics{
private:
	typedef SummaryStatistics self_type;

public:
	SummaryStatistics() :
		total(0),
		total_squared(0),
		n_total(0),
		mean(0),
		standard_deviation(0),
		min(std::numeric_limits<double>::max()),
		max(std::numeric_limits<double>::min())
	{

	}

	bool calculate(void){
		if(this->n_total == 0){
			this->mean = 0;
			this->standard_deviation = 0;
			return false;
		}

		this->mean = this->total / this->n_total;

		if(this->n_total > 1){
			this->standard_deviation = sqrt(this->total_squared/this->n_total - (this->total / this->n_total)*(this->total / this->n_total));
		} else this->standard_deviation = 0;

		return true;
	}

	inline double getSigma(void) const{ return(this->standard_deviation); }
	inline double getSigmaSquared(void) const{ return(this->standard_deviation*this->standard_deviation); }

	template <class T> void operator+=(const T& value){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total       += 1;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	template <class T> void add(const T& value, const double& weight = 1){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total       += weight;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	void reset(void){
		this->total         = 0;
		this->total_squared = 0;
		this->n_total       = 0;
		this->mean          = 0;
		this->standard_deviation = 0;
		this->min = std::numeric_limits<double>::max();
		this->max = std::numeric_limits<double>::min();
	}

	// Accessor functions
	inline double getTotal(void) const{ return(this->total); }
	inline double getTotalSquared(void) const{ return(this->total_squared); }
	inline double getCount(void) const{ return(this->n_total); }
	inline double getMean(void) const{ return(this->mean); }
	inline double getStandardDeviation(void) const{ return(this->standard_deviation); }
	inline double getMin(void) const{ return(this->min); }
	inline double getMax(void) const{ return(this->max); }

	friend std::ostream& operator<<(std::ostream& stream, self_type& self){
		self.calculate();
		stream << self.n_total << "\t" << self.mean << "\t" << self.standard_deviation << "\t" << self.min << "\t" << self.max;
		return(stream);
	}

public:
	double total;
	double total_squared;
	double n_total;
	double mean;
	double standard_deviation;
	double min;
	double max;
};

/**<
 * Standard edge: stores pointers to
 * the start node and the end node
 */
struct Edge{
public:
	typedef Edge               self_type;
	typedef Node               value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;
	typedef SummaryStatistics  score_type;
	typedef tomahawk::io::OutputEntry twk_entry_type;

public:
	Edge(void) :
		hasVisited_(false),
		node_in_(nullptr),
		node_out_(nullptr)
	{

	}

    Edge(const_pointer nodeIn, const_pointer nodeOut) :
    	hasVisited_(false),
		node_in_(nodeIn),
		node_out_(nodeOut)
    {

    }

	Edge(const self_type& other) :
		hasVisited_(other.hasVisited_),
		node_in_(other.node_in_),
		node_out_(other.node_out_)
	{

	}

	~Edge(){}

	inline bool valid(void) const{ return(this->node_in_ != nullptr && this->node_out_ != nullptr && this->isCyclic() == false); }
	inline void isVisited(const bool yes){ this->hasVisited_ = yes; }
	inline bool isCyclic(void) const{ return(this->node_in_ == this->node_out_); }

	/**<
	 * Update the scores of this edge with the provided Tomahawk TWO entry.
	 * This function has its definition in <Node.h>
	 * @param twk_entry Input Tomahawk TWO entry
	 */
	void operator+=(const twk_entry_type& twk_entry);

	void computeStatistics(void){
		this->left_AA_.calculate();
		this->left_AB_.calculate();
		this->left_BA_.calculate();
		this->left_BB_.calculate();
		this->right_AA_.calculate();
		this->right_AB_.calculate();
		this->right_BA_.calculate();
		this->right_BB_.calculate();
	}

	inline const bool hasConnectivity(void) const{
		return(
		       this->left_AA_.n_total + this->left_AB_.n_total + this->left_BA_.n_total + this->left_BB_.n_total +
		       this->right_AA_.n_total + this->right_AB_.n_total + this->right_BA_.n_total + this->right_BB_.n_total
			   );
	}

	score_type& getLargestScore(void){
		score_type* target = &this->left_AA_;
		if(this->left_AB_.n_total > target->n_total) target = &this->left_AB_;
		if(this->left_BA_.n_total > target->n_total) target = &this->left_BA_;
		if(this->left_BB_.n_total > target->n_total) target = &this->left_BB_;
		if(this->right_AB_.n_total > target->n_total) target = &this->right_AB_;
		if(this->right_AB_.n_total > target->n_total) target = &this->right_AB_;
		if(this->right_BA_.n_total > target->n_total) target = &this->right_BA_;
		if(this->right_BB_.n_total > target->n_total) target = &this->right_BB_;
		return(*target);
	}

	inline score_type& getLeftAA(void){ return(this->left_AA_); }
	inline score_type& getLeftAB(void){ return(this->left_AB_); }
	inline score_type& getLeftBA(void){ return(this->left_BA_); }
	inline score_type& getLeftBB(void){ return(this->left_BB_); }
	inline const score_type& getLeftAA(void) const{ return(this->left_AA_); }
	inline const score_type& getLeftAB(void) const{ return(this->left_AB_); }
	inline const score_type& getLeftBA(void) const{ return(this->left_BA_); }
	inline const score_type& getLeftBB(void) const{ return(this->left_BB_); }

	inline score_type& getRightAA(void){ return(this->right_AA_); }
	inline score_type& getRightAB(void){ return(this->right_AB_); }
	inline score_type& getRightBA(void){ return(this->right_BA_); }
	inline score_type& getRightBB(void){ return(this->right_BB_); }
	inline const score_type& getRightAA(void) const{ return(this->right_AA_); }
	inline const score_type& getRightAB(void) const{ return(this->right_AB_); }
	inline const score_type& getRightBA(void) const{ return(this->right_BA_); }
	inline const score_type& getRightBB(void) const{ return(this->right_BB_); }

public:
	bool          hasVisited_;
	const_pointer node_in_;
	const_pointer node_out_;
	score_type    left_AA_, left_AB_, left_BA_, left_BB_;
	score_type    right_AA_, right_AB_, right_BA_, right_BB_;
};

}



#endif /* EDGE_H_ */
