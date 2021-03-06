#ifndef EDGE_CONTAINER_H_
#define EDGE_CONTAINER_H_

namespace javelin{

class EdgeContainer{
public:
	typedef EdgeContainer      self_type;
	typedef Edge               value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef std::ptrdiff_t     difference_type;
    typedef std::size_t        size_type;

public:
    EdgeContainer() :
		n_size_(0),
		n_capacity_(100),
		edges_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{

	}

    EdgeContainer(const size_t capacity) :
		n_size_(0),
		n_capacity_(capacity),
		edges_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{

	}

    EdgeContainer(const self_type& other) :
		n_size_(other.n_size_),
		n_capacity_(other.n_capacity_),
		edges_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{
		for(size_t i = 0; i < this->size(); ++i)
			new( &this->edges_[i] ) value_type( other[i] );
	}

    ~EdgeContainer(){
		for(std::size_t i = 0; i < this->size(); ++i)
			((this->edges_ + i)->~Edge)();

		::operator delete[](static_cast<void*>(this->edges_));
	}

    class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->edges_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->edges_[position]); }
	inline reference operator[](const size_type& position){ return(this->edges_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->edges_[position]); }
	inline pointer data(void){ return(this->edges_); }
	inline const_pointer data(void) const{ return(this->edges_); }
	inline reference front(void){ return(this->edges_[0]); }
	inline const_reference front(void) const{ return(this->edges_[0]); }
	inline reference back(void){ return(this->edges_[this->n_size_ - 1]); }
	inline const_reference back(void) const{ return(this->edges_[this->n_size_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_size_ == 0); }
	inline const size_type& size(void) const{ return(this->n_size_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->edges_[0]); }
	inline iterator end()  { return iterator(&this->edges_[this->n_size_]); }
	inline const_iterator begin()  const{ return const_iterator(&this->edges_[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->edges_[this->n_size_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->edges_[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->edges_[this->n_size_]); }

	// Overload
	inline self_type& operator+=(const_reference edge){
		if(this->size() + 1 == this->capacity()){
			std::cerr << "is full -> resize: " << this->size()+1 << "/" << this->capacity() << std::endl;
			this->resize();
			return(*this);
		}
		new( &this->edges_[this->n_size_] ) value_type( edge );
		++this->n_size_;
		return(*this);
	}

	void resize(const size_t new_size){
		if(new_size <= this->capacity()){
			// Delete data from `new_size` until `this->size()`
			for(size_t i = new_size; i < this->size(); ++i)
				((this->edges_ + i)->~Edge)();

			this->n_size_ = new_size;
			return;
		} else {
			this->n_capacity_ = new_size;
		}

		pointer old = this->edges_;
		this->edges_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
		for(size_t i = 0; i < this->size(); ++i){
			new( &this->edges_[i] ) value_type( old[i] );
			((old + i)->~Edge)();
		}
		::operator delete[](static_cast<void*>(old));
	}
	inline void resize(void){ this->resize(this->capacity()*2); }

public:
    size_type  n_size_;
    size_type  n_capacity_;
    pointer    edges_;
};

}


#endif /* EDGE_CONTAINER_H_ */
