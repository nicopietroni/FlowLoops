#ifndef _UTILITY_FUNCTIONS_

#define _UTILITY_FUNCTIONS_

#ifndef GX_DEBUG_ASSERT
#ifdef NDEBUG  // if not in debug mode
    #define gx_debug_assert(cond) ((void)0)
#else  // debug mode
    #define gx_debug_assert(cond) assert(cond)
#endif
#endif

#ifndef GX_ASSERT_NOT_REACHED
#ifdef NDEBUG  // if not in debug mode
    #define gx_assert_not_reached() ((void)0)  // do nothing in release builds
#else  // debug mode
    #define gx_assert_not_reached() assert(false)  // assert false, indicating this code should not be reached
#endif
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::Vector3d;
using Eigen::Vector4d;

#include <vector>
#include <set>

#include <iostream>

#include "cs_logging.h"

	template <class T>
	inline T sign(const T& val){
		if(val >= T(0)){
			return T(1);
		}
		else 
			return T(-1);
	}

	template <class T>
	inline T mod(const T& val, const T& base){
		if(base <= T(0))
		{return 0;}
		else{

			T ret = val%base;

			if (ret < T(0))
			{
				return ret + base;
			}
			else
			{
				return ret;
			}
		}
	}

	template <class T>
	class MyPair {
	public:
		T id[2];
		MyPair(){
			id[0] = id[1] = T(-1);
		}
		MyPair(T a, T b){
			id[0] = a;
			id[1] = b;
		}
	};


	
	//------------------------------ Small set ----------------------------------------------

	  /**
     * Similar to std::set, but with fixed maximum size
     * (and no dynamic memory allocation). Used to store
     * vertices equations (represented as plane indices 
     * triplets).
     */
    template <class T, int DIM> class small_set {
    public:
        typedef small_set<T,DIM> thisclass ;
        typedef T* iterator ;
        typedef const T* const_iterator ;
        typedef T& reference ;
        typedef T value_type ;
        
        small_set() : end_(data_) { }

        small_set(const thisclass& rhs) {
            copy(rhs) ;
        }
        
        thisclass& operator=(const thisclass& rhs) {
            copy(rhs) ; return *this ;
        }

        template <int DIM2> small_set(const small_set<T,DIM2>& rhs) {
            copy(rhs) ;
        }
        
        template <int DIM2> thisclass& operator=(const small_set<T,DIM2>& rhs) {
            copy(rhs) ; return *this ;
        }

        unsigned int size() const { return (unsigned int)(end_ - data_) ; }
        unsigned int capacity() const { return (unsigned int)DIM ; }

        iterator begin() { return data_ ; }
        iterator end() { return end_ ; }
        iterator end_of_storage() { return data_ + DIM ; }

        iterator* end_ptr() { return &end_ ; }

        const_iterator begin() const { return data_ ; }
        const_iterator end() const { return end_ ; }
        const_iterator end_of_storage() const { return data_ + DIM ; }

        iterator insert(const T& x) { 
            return insert(x, find_i(x)) ; 
        }

        iterator insert(const T& x, iterator where) {
            if(where == end()) { *where = x ; grow() ; return where ; }
            if(*where == x) { return where ; }
            grow() ;
            if(where == end() - 1) {
				Geex::cs_err() << ">>>>>>>>Error out of set bound<<<<<<<<" << std::endl; //hpan
                *where = x ; return where ; 
            }
            for(iterator i = end()-1; i != where; i--) {
                gx_debug_assert(i != begin()) ;
                *i = *(i-1) ;
            }
            *where = x ;
#ifdef GX_DEBUG
            for(iterator i = begin(); i != end() - 1; i++) {
                gx_debug_assert(*i < *(i+1)) ;
            }
#endif
            return where ;
        }

		//<<<<<
		void insert(const thisclass& rhs) {
			for (const_iterator i = rhs.begin(); i != rhs.end(); ++i)
			{
				this->insert(*i, find_i(*i));
			}
		}

		inline bool operator<(const thisclass& rhs) const{
			int com_size = (size() < rhs.size()) ? size() : rhs.size();
			for (int i = 0; i < com_size; i++)
			{
				if ((*this)[i] < rhs[i])
				{
					return true;
				}
				else
					if ((*this)[i] > rhs[i])
					{
						return false;
					}
			}
			if (size() < rhs.size())
				return true;
			return false;
		}

		inline bool operator==(const thisclass& rhs) const{
			if (size() != rhs.size())
			{
				return false;
			}
			
			for (int i = 0; i < size(); i++)
			{
				if ((*this)[i] != rhs[i])
				{
					return false;
				}
			}

			return true;
		}

		//>>>>>

        void clear() { end_ = data_ ; }

        iterator find(const T& x) {
            iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        const_iterator find(const T& x) const {
            const_iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        void push_back(const T& x) {
#ifdef GX_DEBUG
            for(iterator i = begin(); i != end(); i++) {
                gx_debug_assert(*i < x) ;
            }
#endif
            *end_ = x ;
            grow() ;
        }

        void print(std::ostream& out) const {
            out << "[ " ;
            for(const_iterator it=begin(); it!=end(); it++) {
                out << *it << " " ;
            }
            out << "]" ;
        }

        T& operator[](int i) {
            gx_debug_assert(i >= 0) ;
            gx_debug_assert(begin() + i < end()) ;
            return begin()[i] ;
        }

        const T& operator[](int i) const {
            gx_debug_assert(i >= 0) ;
            gx_debug_assert(begin() + i < end()) ;
            return begin()[i] ;
        }

    protected:

        void grow() {
            gx_debug_assert(end() != end_of_storage()) ;
            end_++ ;
        }


        template <int DIM2> void copy(const small_set<T,DIM2>& rhs) {
            end_ = data_ ;
            for(typename small_set<T,DIM2>::const_iterator it=rhs.begin(); it!=rhs.end(); it++) {
                push_back(*it) ;
            }
        }


        // Note: maybe we should start from end() instead of begin()
        // since negative indices are inserted first.
        iterator find_i(const T& x) {
            iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

        const_iterator find_i(const T& x) const {
            const_iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

    protected:
        T data_[DIM] ;
        iterator end_ ;
    } ;

    template <class T, int DIM> inline std::ostream& operator<<(
        std::ostream& out,
        const small_set<T, DIM>& S
    ) {
        S.print(out) ; return out ;
    }


    template <class T, int DIM> inline std::ostream& operator<<(
        std::ostream& out,
        const std::set<T>& S
    ) {
        out << "[ " ;
        for(typename std::set<T>::const_iterator it = S.begin(); it != S.end(); it++) {
            out << *it << " " ;
        }
        out << "]" ;
        return out ;
    }
    
	//--------------------------------------------------------------------------

#endif