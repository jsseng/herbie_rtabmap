/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
 * Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
 *
 * THE BSD LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#ifndef RTABMAP_FLANN_KDTREE_INDEX_H_
#define RTABMAP_FLANN_KDTREE_INDEX_H_

#include <algorithm>
#include <map>
#include <cassert>
#include <cstring>
#include <stdarg.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <sys/mman.h>

#include "rtflann/general.h"
#include "rtflann/algorithms/nn_index.h"
#include "rtflann/util/dynamic_bitset.h"
#include "rtflann/util/matrix.h"
#include "rtflann/util/result_set.h"
#include "rtflann/util/heap.h"
#include "rtflann/util/allocator.h"
#include "rtflann/util/random.h"
#include "rtflann/util/saving.h"


using std::cout;
using std::endl;
using std::vector;

namespace rtflann
{

struct KDTreeIndexParams : public IndexParams
{
    KDTreeIndexParams(int trees = 4)
    {
        (*this)["algorithm"] = FLANN_INDEX_KDTREE;
        (*this)["trees"] = trees;
    }
};


/**
 * Randomized kd-tree index
 *
 * Contains the k-d trees and other information for indexing a set of points
 * for nearest-neighbor matching.
 */
template <typename Distance>
class KDTreeIndex : public NNIndex<Distance>
{
public:
    typedef typename Distance::ElementType ElementType;
    typedef typename Distance::ResultType DistanceType;

    typedef NNIndex<Distance> BaseClass;

    typedef bool needs_kdtree_distance;

private:
	 /*--------------------- Internal Data Structures --------------------------*/
	struct Node
	{
		/**
		 * Dimension used for subdivision.
		 */
		int divfeat;
		/**
		 * The values used for subdivision.
		 */
		DistanceType divval;
		/**
		 * Point data
		 */
		ElementType* point;
        int point_num;
		/**
		* The child nodes.
		*/
		Node* child1, *child2;
        int fix;
		Node(){
			child1 = NULL;
			child2 = NULL;
		}
		~Node() {
			if (child1 != NULL) { child1->~Node(); child1 = NULL; }

			if (child2 != NULL) { child2->~Node(); child2 = NULL; }
		}

	private:
		template<typename Archive>
		void serialize(Archive& ar)
		{
			typedef KDTreeIndex<Distance> Index;
			Index* obj = static_cast<Index*>(ar.getObject());

			ar & divfeat;
			ar & divval;

			bool leaf_node = false;
			if (Archive::is_saving::value) {
				leaf_node = ((child1==NULL) && (child2==NULL));
			}
			ar & leaf_node;

			if (leaf_node) {
				if (Archive::is_loading::value) {
					point = obj->points_[divfeat];
				}
			}

			if (!leaf_node) {
				if (Archive::is_loading::value) {
					child1 = new(obj->pool_) Node();
					child2 = new(obj->pool_) Node();
				}
				ar & *child1;
				ar & *child2;
			}
		}
		friend struct serialization::access;
	};

	typedef Node* NodePtr;
	typedef BranchStruct<NodePtr, DistanceType> BranchSt;
	typedef BranchSt* Branch;

public:

    /**
     * KDTree constructor
     *
     * Params:
     *          inputData = dataset with the input features
     *          params = parameters passed to the kdtree algorithm
     */
    KDTreeIndex(const IndexParams& params = KDTreeIndexParams(), Distance d = Distance() ) :
    	BaseClass(params, d), mean_(NULL), var_(NULL)
    {
        trees_ = get_param(index_params_,"trees",4);
        load_cached = 0;
    }


    /**
     * KDTree constructor
     *
     * Params:
     *          inputData = dataset with the input features
     *          params = parameters passed to the kdtree algorithm
     */
    KDTreeIndex(const Matrix<ElementType>& dataset, const IndexParams& params = KDTreeIndexParams(),
                Distance d = Distance() ) : BaseClass(params,d ), mean_(NULL), var_(NULL)
    {
        trees_ = get_param(index_params_,"trees",4);
        load_cached = 0;
        this->cached = 0;

        setDataset(dataset);
    }

    KDTreeIndex(const KDTreeIndex& other) : BaseClass(other),
    		trees_(other.trees_)
    {
        tree_roots_.resize(other.tree_roots_.size());
        for (size_t i=0;i<tree_roots_.size();++i) {
        	copyTree(tree_roots_[i], other.tree_roots_[i]);
        }
        load_cached = 0;
    }

    KDTreeIndex& operator=(KDTreeIndex other)
    {
    	this->swap(other);
    	return *this;
    }

    /**
     * Standard destructor
     */
    virtual ~KDTreeIndex()
    {
    	freeIndex();
    }

    BaseClass* clone() const
    {
    	return new KDTreeIndex(*this);
    }

    using BaseClass::buildIndex;
    
    void addPoints(const Matrix<ElementType>& points, float rebuild_threshold = 2)
    {
        assert(points.cols==veclen_);

        size_t old_size = size_;
        //std::cout << "---(addPoints() - kdtree_index.h)  size: " << old_size << "----" << std::endl;
        extendDataset(points);
        
        if (rebuild_threshold>1 && size_at_build_*rebuild_threshold<size_) {
            buildIndex();
        }
        else {
            for (size_t i=old_size;i<size_;++i) {
                for (int j = 0; j < trees_; j++) {
                    addPointToTree(tree_roots_[j], i);
                }
            }
        }        
    }

    flann_algorithm_t getType() const
    {
        return FLANN_INDEX_KDTREE;
    }


    template<typename Archive>
    void serialize(Archive& ar)
    {
    	ar.setObject(this);

    	ar & *static_cast<NNIndex<Distance>*>(this);

    	ar & trees_;

    	if (Archive::is_loading::value) {
    		tree_roots_.resize(trees_);
    	}
    	for (size_t i=0;i<tree_roots_.size();++i) {
    		if (Archive::is_loading::value) {
    			tree_roots_[i] = new(pool_) Node();
    		}
    		ar & *tree_roots_[i];
    	}

    	if (Archive::is_loading::value) {
            index_params_["algorithm"] = getType();
            index_params_["trees"] = trees_;
    	}
    }


    void saveIndex(FILE* stream)
    {
    	serialization::SaveArchive sa(stream);
    	sa & *this;
    }


    void loadIndex(FILE* stream)
    {
    	freeIndex();
    	serialization::LoadArchive la(stream);
    	la & *this;
    }

    /**
     * Computes the inde memory usage
     * Returns: memory used by the index
     */
    int usedMemory() const
    {
        return int(pool_.usedMemory+pool_.wastedMemory+size_*sizeof(int));  // pool memory and vind array memory
    }

    /**
     * Find set of nearest neighbors to vec. Their indices are stored inside
     * the result object.
     *
     * Params:
     *     result = the result object in which the indices of the nearest-neighbors are stored
     *     vec = the vector for which to search the nearest neighbors
     *     maxCheck = the maximum number of restarts (in a best-bin-first manner)
     */
    void findNeighbors(ResultSet<DistanceType>& result, const ElementType* vec, const SearchParams& searchParams) const
    {
        int maxChecks = searchParams.checks;
        float epsError = 1+searchParams.eps;

        if (maxChecks==FLANN_CHECKS_UNLIMITED) {
        	if (removed_) {
        		getExactNeighbors<true>(result, vec, epsError);
        	}
        	else {
        		getExactNeighbors<false>(result, vec, epsError);
        	}
        }
        else {
        	if (removed_) {
        		getNeighbors<true>(result, vec, maxChecks, epsError);
        	}
        	else {
        		getNeighbors<false>(result, vec, maxChecks, epsError);
        	}
        }
    }

#ifdef FLANN_KDTREE_MEM_OPT

    /**
	 * Find set of nearest neighbors to vec. Their indices are stored inside
	 * the result object.
	 *
	 * Params:
	 *     result = the result object in which the indices of the nearest-neighbors are stored
	 *     vec = the vector for which to search the nearest neighbors
	 *     maxCheck = the maximum number of restarts (in a best-bin-first manner)
	 */
	void findNeighbors(ResultSet<DistanceType>& result, const ElementType* vec, const SearchParams& searchParams, Heap<BranchSt>* heap) const
	{
		int maxChecks = searchParams.checks;
		float epsError = 1+searchParams.eps;

		if (maxChecks==FLANN_CHECKS_UNLIMITED) {
			if (removed_) {
				getExactNeighbors<true>(result, vec, epsError);
			}
			else {
				getExactNeighbors<false>(result, vec, epsError);
			}
		}
		else {
			if (removed_) {
				getNeighbors<true>(result, vec, maxChecks, epsError, heap);
			}
			else {
				getNeighbors<false>(result, vec, maxChecks, epsError, heap);
			}
		}
	}

	/**
	 * @brief Perform k-nearest neighbor search
	 * @param[in] queries The query points for which to find the nearest neighbors
	 * @param[out] indices The indices of the nearest neighbors found
	 * @param[out] dists Distances to the nearest neighbors found
	 * @param[in] knn Number of nearest neighbors to return
	 * @param[in] params Search parameters
	 */
	virtual int knnSearch(const Matrix<ElementType>& queries,
			Matrix<size_t>& indices,
			Matrix<DistanceType>& dists,
			size_t knn,
			const SearchParams& params) const
	{
		assert(queries.cols == veclen());
		assert(indices.rows >= queries.rows);
		assert(dists.rows >= queries.rows);
		assert(indices.cols >= knn);
		assert(dists.cols >= knn);
		bool use_heap;

		if (params.use_heap==FLANN_Undefined) {
			use_heap = (knn>KNN_HEAP_THRESHOLD)?true:false;
		}
		else {
			use_heap = (params.use_heap==FLANN_True)?true:false;
		}
		int count = 0;

		Heap<BranchSt>* heap = new Heap<BranchSt>((int)size_);

		if (use_heap) {
	//#pragma omp parallel num_threads(params.cores)
			{
				KNNResultSet2<DistanceType> resultSet(knn);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					size_t n = std::min(resultSet.size(), knn);
					resultSet.copy(indices[i], dists[i], n, params.sorted);
					indices_to_ids(indices[i], indices[i], n);
					count += n;
				}
			}
		}
		else {
			std::vector<double> times(queries.rows);
	//#pragma omp parallel num_threads(params.cores)
			{
				KNNSimpleResultSet<DistanceType> resultSet(knn);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					size_t n = std::min(resultSet.size(), knn);
					resultSet.copy(indices[i], dists[i], n, params.sorted);
					indices_to_ids(indices[i], indices[i], n);
					count += n;
				}
			}
			std::sort(times.begin(), times.end());
		}
		delete heap;
		return count;
	}


	/**
	 * @brief Perform k-nearest neighbor search
	 * @param[in] queries The query points for which to find the nearest neighbors
	 * @param[out] indices The indices of the nearest neighbors found
	 * @param[out] dists Distances to the nearest neighbors found
	 * @param[in] knn Number of nearest neighbors to return
	 * @param[in] params Search parameters
	 */
	virtual int knnSearch(const Matrix<ElementType>& queries,
					std::vector< std::vector<size_t> >& indices,
					std::vector<std::vector<DistanceType> >& dists,
					size_t knn,
					const SearchParams& params) const
	{
		assert(queries.cols == veclen());
		bool use_heap;
		if (params.use_heap==FLANN_Undefined) {
			use_heap = (knn>KNN_HEAP_THRESHOLD)?true:false;
		}
		else {
			use_heap = (params.use_heap==FLANN_True)?true:false;
		}

		if (indices.size() < queries.rows ) indices.resize(queries.rows);
		if (dists.size() < queries.rows ) dists.resize(queries.rows);

		Heap<BranchSt>* heap = new Heap<BranchSt>((int)size_);

		int count = 0;
		if (use_heap) {
	//#pragma omp parallel num_threads(params.cores)
			{
				KNNResultSet2<DistanceType> resultSet(knn);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					size_t n = std::min(resultSet.size(), knn);
					indices[i].resize(n);
					dists[i].resize(n);
					if (n>0) {
						resultSet.copy(&indices[i][0], &dists[i][0], n, params.sorted);
						indices_to_ids(&indices[i][0], &indices[i][0], n);
					}
					count += n;
				}
			}
		}
		else {
	//#pragma omp parallel num_threads(params.cores)
			{
				KNNSimpleResultSet<DistanceType> resultSet(knn);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					size_t n = std::min(resultSet.size(), knn);
					indices[i].resize(n);
					dists[i].resize(n);
					if (n>0) {
						resultSet.copy(&indices[i][0], &dists[i][0], n, params.sorted);
						indices_to_ids(&indices[i][0], &indices[i][0], n);
					}
					count += n;
				}
			}
		}
		delete heap;

		return count;
	}

	/**
	 * @brief Perform radius search
	 * @param[in] query The query point
	 * @param[out] indices The indices of the neighbors found within the given radius
	 * @param[out] dists The distances to the nearest neighbors found
	 * @param[in] radius The radius used for search
	 * @param[in] params Search parameters
	 * @return Number of neighbors found
	 */
	virtual int radiusSearch(const Matrix<ElementType>& queries,
			Matrix<size_t>& indices,
			Matrix<DistanceType>& dists,
			float radius,
			const SearchParams& params) const
	{
		assert(queries.cols == veclen());
		int count = 0;
		size_t num_neighbors = std::min(indices.cols, dists.cols);
		int max_neighbors = params.max_neighbors;
		if (max_neighbors<0) max_neighbors = num_neighbors;
		else max_neighbors = std::min(max_neighbors,(int)num_neighbors);

		Heap<BranchSt>* heap = new Heap<BranchSt>((int)size_);

		if (max_neighbors==0) {
	//#pragma omp parallel num_threads(params.cores)
			{
				CountRadiusResultSet<DistanceType> resultSet(radius);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					count += resultSet.size();
				}
			}
		}
		else {
			// explicitly indicated to use unbounded radius result set
			// and we know there'll be enough room for resulting indices and dists
			if (params.max_neighbors<0 && (num_neighbors>=this->size())) {
	//#pragma omp parallel num_threads(params.cores)
				{
					RadiusResultSet<DistanceType> resultSet(radius);
	//#pragma omp for schedule(static) reduction(+:count)
					for (int i = 0; i < (int)queries.rows; i++) {
						resultSet.clear();
						findNeighbors(resultSet, queries[i], params, heap);
						size_t n = resultSet.size();
						count += n;
						if (n>num_neighbors) n = num_neighbors;
						resultSet.copy(indices[i], dists[i], n, params.sorted);

						// mark the next element in the output buffers as unused
						if (n<indices.cols) indices[i][n] = size_t(-1);
						if (n<dists.cols) dists[i][n] = std::numeric_limits<DistanceType>::infinity();
						indices_to_ids(indices[i], indices[i], n);
					}
				}
			}
			else {
				// number of neighbors limited to max_neighbors
	//#pragma omp parallel num_threads(params.cores)
				{
					KNNRadiusResultSet<DistanceType> resultSet(radius, max_neighbors);
	//#pragma omp for schedule(static) reduction(+:count)
					for (int i = 0; i < (int)queries.rows; i++) {
						resultSet.clear();
						findNeighbors(resultSet, queries[i], params, heap);
						size_t n = resultSet.size();
						count += n;
						if ((int)n>max_neighbors) n = max_neighbors;
						resultSet.copy(indices[i], dists[i], n, params.sorted);

						// mark the next element in the output buffers as unused
						if (n<indices.cols) indices[i][n] = size_t(-1);
						if (n<dists.cols) dists[i][n] = std::numeric_limits<DistanceType>::infinity();
						indices_to_ids(indices[i], indices[i], n);
					}
				}
			}
		}
		delete heap;
		return count;
	}

	/**
	 * @brief Perform radius search
	 * @param[in] query The query point
	 * @param[out] indices The indices of the neighbors found within the given radius
	 * @param[out] dists The distances to the nearest neighbors found
	 * @param[in] radius The radius used for search
	 * @param[in] params Search parameters
	 * @return Number of neighbors found
	 */
	virtual int radiusSearch(const Matrix<ElementType>& queries,
			std::vector< std::vector<size_t> >& indices,
			std::vector<std::vector<DistanceType> >& dists,
			float radius,
			const SearchParams& params) const
	{
		assert(queries.cols == veclen());
		int count = 0;

		Heap<BranchSt>* heap = new Heap<BranchSt>((int)size_);

		// just count neighbors
		if (params.max_neighbors==0) {
	//#pragma omp parallel num_threads(params.cores)
			{
				CountRadiusResultSet<DistanceType> resultSet(radius);
	//#pragma omp for schedule(static) reduction(+:count)
				for (int i = 0; i < (int)queries.rows; i++) {
					resultSet.clear();
					findNeighbors(resultSet, queries[i], params, heap);
					count += resultSet.size();
				}
			}
		}
		else {
			if (indices.size() < queries.rows ) indices.resize(queries.rows);
			if (dists.size() < queries.rows ) dists.resize(queries.rows);

			if (params.max_neighbors<0) {
				// search for all neighbors
	//#pragma omp parallel num_threads(params.cores)
				{
					RadiusResultSet<DistanceType> resultSet(radius);
	//#pragma omp for schedule(static) reduction(+:count)
					for (int i = 0; i < (int)queries.rows; i++) {
						resultSet.clear();
						findNeighbors(resultSet, queries[i], params, heap);
						size_t n = resultSet.size();
						count += n;
						indices[i].resize(n);
						dists[i].resize(n);
						if (n > 0) {
							resultSet.copy(&indices[i][0], &dists[i][0], n, params.sorted);
							indices_to_ids(&indices[i][0], &indices[i][0], n);
						}
					}
				}
			}
			else {
				// number of neighbors limited to max_neighbors
	//#pragma omp parallel num_threads(params.cores)
				{
					KNNRadiusResultSet<DistanceType> resultSet(radius, params.max_neighbors);
	//#pragma omp for schedule(static) reduction(+:count)
					for (int i = 0; i < (int)queries.rows; i++) {
						resultSet.clear();
						findNeighbors(resultSet, queries[i], params, heap);
						size_t n = resultSet.size();
						count += n;
						if ((int)n>params.max_neighbors) n = params.max_neighbors;
						indices[i].resize(n);
						dists[i].resize(n);
						if (n > 0) {
							resultSet.copy(&indices[i][0], &dists[i][0], n, params.sorted);
							indices_to_ids(&indices[i][0], &indices[i][0], n);
						}
					}
				}
			}
		}
		delete heap;
		return count;
	}
#endif

protected:

    /**
     * Builds the index
     */
    void buildIndexImpl()
    {
        if (this->cached==0) {
            std::cout << "---(buildIndexImpl() - kdtree_index.h)---" << std::endl;
            // Create a permutable array of indices to the input vectors.
            std::vector<int> ind(size_);
            for (size_t i = 0; i < size_; ++i)
            {
                ind[i] = int(i);
            }

            mean_ = new DistanceType[veclen_];
            var_ = new DistanceType[veclen_];

            tree_roots_.resize(trees_);
            /* Construct the randomized trees. */
            for (int i = 0; i < trees_; i++)
            {
                /* Randomize the order of vectors to allow for unbiased sampling. */
                std::random_shuffle(ind.begin(), ind.end());
                tree_roots_[i] = divideTree(&ind[0], int(size_));
            }
            delete[] mean_;
            delete[] var_;
        }
    }

    void freeIndex()
    {
        //if (load_cached == 0) {
            for (size_t i = 0; i < tree_roots_.size(); ++i)
            {
                // using placement new, so call destructor explicitly
                if (tree_roots_[i] != NULL)
                    tree_roots_[i]->~Node();
            }
            pool_.free();
        //}
    }


private:

    void copyTree(NodePtr& dst, const NodePtr& src)
    {
    	dst = new(pool_) Node();
    	dst->divfeat = src->divfeat;
    	dst->divval = src->divval;
    	if (src->child1==NULL && src->child2==NULL) {
    		dst->point = points_[dst->divfeat];
            dst->point_num = dst->divfeat;
    		dst->child1 = NULL;
    		dst->child2 = NULL;
    	}
    	else {
    		copyTree(dst->child1, src->child1);
    		copyTree(dst->child2, src->child2);
    	}
    }

    /**
     * Create a tree node that subdivides the list of vecs from vind[first]
     * to vind[last].  The routine is called recursively on each sublist.
     * Place a pointer to this new tree node in the location pTree.
     *
     * Params: pTree = the new node to create
     *                  first = index of the first vector
     *                  last = index of the last vector
     */
    NodePtr divideTree(int* ind, int count)
    {
        NodePtr node = new(pool_) Node(); // allocate memory

        /* If too few exemplars remain, then make this a leaf node. */
        if (count == 1) {
            node->child1 = node->child2 = NULL;    /* Mark as leaf node. */
            node->divfeat = *ind;    /* Store index of this vec. */
            node->point = points_[*ind];
            node->point_num = *ind;
        }
        else {
            int idx;
            int cutfeat;
            DistanceType cutval;
            meanSplit(ind, count, idx, cutfeat, cutval);

            node->divfeat = cutfeat;
            node->divval = cutval;
            node->child1 = divideTree(ind, idx);
            node->child2 = divideTree(ind+idx, count-idx);
        }

        return node;
    }


    /**
     * Choose which feature to use in order to subdivide this set of vectors.
     * Make a random choice among those with the highest variance, and use
     * its variance as the threshold value.
     */
    void meanSplit(int* ind, int count, int& index, int& cutfeat, DistanceType& cutval)
    {
        memset(mean_,0,veclen_*sizeof(DistanceType));
        memset(var_,0,veclen_*sizeof(DistanceType));

        /* Compute mean values.  Only the first SAMPLE_MEAN values need to be
            sampled to get a good estimate.
         */
        int cnt = std::min((int)SAMPLE_MEAN+1, count);
        for (int j = 0; j < cnt; ++j) {
            ElementType* v = points_[ind[j]];
            for (size_t k=0; k<veclen_; ++k) {
                mean_[k] += v[k];
            }
        }
        DistanceType div_factor = DistanceType(1)/cnt;
        for (size_t k=0; k<veclen_; ++k) {
            mean_[k] *= div_factor;
        }

        /* Compute variances (no need to divide by count). */
        for (int j = 0; j < cnt; ++j) {
            ElementType* v = points_[ind[j]];
            for (size_t k=0; k<veclen_; ++k) {
                DistanceType dist = v[k] - mean_[k];
                var_[k] += dist * dist;
            }
        }
        /* Select one of the highest variance indices at random. */
        cutfeat = selectDivision(var_);
        cutval = mean_[cutfeat];

        int lim1, lim2;
        planeSplit(ind, count, cutfeat, cutval, lim1, lim2);

        if (lim1>count/2) index = lim1;
        else if (lim2<count/2) index = lim2;
        else index = count/2;

        /* If either list is empty, it means that all remaining features
         * are identical. Split in the middle to maintain a balanced tree.
         */
        if ((lim1==count)||(lim2==0)) index = count/2;
    }


    /**
     * Select the top RAND_DIM largest values from v and return the index of
     * one of these selected at random.
     */
    int selectDivision(DistanceType* v)
    {
        int num = 0;
        size_t topind[RAND_DIM];

        /* Create a list of the indices of the top RAND_DIM values. */
        for (size_t i = 0; i < veclen_; ++i) {
            if ((num < RAND_DIM)||(v[i] > v[topind[num-1]])) {
                /* Put this element at end of topind. */
                if (num < RAND_DIM) {
                    topind[num++] = i;            /* Add to list. */
                }
                else {
                    topind[num-1] = i;         /* Replace last element. */
                }
                /* Bubble end value down to right location by repeated swapping. */
                int j = num - 1;
                while (j > 0  &&  v[topind[j]] > v[topind[j-1]]) {
                    std::swap(topind[j], topind[j-1]);
                    --j;
                }
            }
        }
        /* Select a random integer in range [0,num-1], and return that index. */
        int rnd = rand_int(num);
        return (int)topind[rnd];
    }


    /**
     *  Subdivide the list of points by a plane perpendicular on axe corresponding
     *  to the 'cutfeat' dimension at 'cutval' position.
     *
     *  On return:
     *  dataset[ind[0..lim1-1]][cutfeat]<cutval
     *  dataset[ind[lim1..lim2-1]][cutfeat]==cutval
     *  dataset[ind[lim2..count]][cutfeat]>cutval
     */
    void planeSplit(int* ind, int count, int cutfeat, DistanceType cutval, int& lim1, int& lim2)
    {
        /* Move vector indices for left subtree to front of list. */
        int left = 0;
        int right = count-1;
        for (;; ) {
            while (left<=right && points_[ind[left]][cutfeat]<cutval) ++left;
            while (left<=right && points_[ind[right]][cutfeat]>=cutval) --right;
            if (left>right) break;
            std::swap(ind[left], ind[right]); ++left; --right;
        }
        lim1 = left;
        right = count-1;
        for (;; ) {
            while (left<=right && points_[ind[left]][cutfeat]<=cutval) ++left;
            while (left<=right && points_[ind[right]][cutfeat]>cutval) --right;
            if (left>right) break;
            std::swap(ind[left], ind[right]); ++left; --right;
        }
        lim2 = left;
    }

    /**
     * Performs an exact nearest neighbor search. The exact search performs a full
     * traversal of the tree.
     */
    template<bool with_removed>
    void getExactNeighbors(ResultSet<DistanceType>& result, const ElementType* vec, float epsError) const
    {
        //		checkID -= 1;  /* Set a different unique ID for each search. */

        if (trees_ > 1) {
            fprintf(stderr,"It doesn't make any sense to use more than one tree for exact search");
        }
        if (trees_>0) {
            searchLevelExact<with_removed>(result, vec, tree_roots_[0], 0.0, epsError);
        }
    }

    /**
     * Performs the approximate nearest-neighbor search. The search is approximate
     * because the tree traversal is abandoned after a given number of descends in
     * the tree.
     */
    template<bool with_removed>
    void getNeighbors(ResultSet<DistanceType>& result, const ElementType* vec, int maxCheck, float epsError) const
    {
        int i;
        BranchSt branch;

        int checkCount = 0;
        Heap<BranchSt>* heap = new Heap<BranchSt>((int)size_);
        DynamicBitset checked(size_);

        /* Search once through each tree down to root. */
        for (i = 0; i < trees_; ++i) {
            searchLevel<with_removed>(result, vec, tree_roots_[i], 0, checkCount, maxCheck, epsError, heap, checked);
        }

        /* Keep searching other branches from heap until finished. */
        while ( heap->popMin(branch) && (checkCount < maxCheck || !result.full() )) {
            searchLevel<with_removed>(result, vec, branch.node, branch.mindist, checkCount, maxCheck, epsError, heap, checked);
        }

        delete heap;

    }

#ifdef FLANN_KDTREE_MEM_OPT
    /**
	 * Performs the approximate nearest-neighbor search. The search is approximate
	 * because the tree traversal is abandoned after a given number of descends in
	 * the tree.
	 */
	template<bool with_removed>
	void getNeighbors(ResultSet<DistanceType>& result, const ElementType* vec, int maxCheck, float epsError, Heap<BranchSt>* heap) const
	{
		int i;
		BranchSt branch;

		int checkCount = 0;
		DynamicBitset checked(size_);
		heap->clear();

		/* Search once through each tree down to root. */
		for (i = 0; i < trees_; ++i) {
			searchLevel<with_removed>(result, vec, tree_roots_[i], 0, checkCount, maxCheck, epsError, heap, checked);
		}

		/* Keep searching other branches from heap until finished. */
		while ( heap->popMin(branch) && (checkCount < maxCheck || !result.full() )) {
			searchLevel<with_removed>(result, vec, branch.node, branch.mindist, checkCount, maxCheck, epsError, heap, checked);
		}
	}
#endif


    /**
     *  Search starting from a given node of the tree.  Based on any mismatches at
     *  higher levels, all exemplars below this level must have a distance of
     *  at least "mindistsq".
     */
    template<bool with_removed>
    void searchLevel(ResultSet<DistanceType>& result_set, const ElementType* vec, NodePtr node, DistanceType mindist, int& checkCount, int maxCheck,
                     float epsError, Heap<BranchSt>* heap, DynamicBitset& checked) const
    {
        if (result_set.worstDist()<mindist) {
            //			printf("Ignoring branch, too far\n");
            return;
        }

        /* If this is a leaf node, then do check and return. */
        if ((node->child1 == NULL)&&(node->child2 == NULL)) {
            int index = node->divfeat;
            if (with_removed) {
            	if (removed_points_.test(index)) return;
            }
            /*  Do not check same node more than once when searching multiple trees. */
            if ( checked.test(index) || ((checkCount>=maxCheck)&& result_set.full()) ) return;
            checked.set(index);
            checkCount++;

            DistanceType dist = distance_(node->point, vec, veclen_);
            result_set.addPoint(dist,index);
            return;
        }

        /* Which child branch should be taken first? */
        ElementType val = vec[node->divfeat];
        DistanceType diff = val - node->divval;
        NodePtr bestChild = (diff < 0) ? node->child1 : node->child2;
        NodePtr otherChild = (diff < 0) ? node->child2 : node->child1;

        /* Create a branch record for the branch not taken.  Add distance
            of this feature boundary (we don't attempt to correct for any
            use of this feature in a parent node, which is unlikely to
            happen and would have only a small effect).  Don't bother
            adding more branches to heap after halfway point, as cost of
            adding exceeds their value.
         */

        DistanceType new_distsq = mindist + distance_.accum_dist(val, node->divval, node->divfeat);
        //		if (2 * checkCount < maxCheck  ||  !result.full()) {
        if ((new_distsq*epsError < result_set.worstDist())||  !result_set.full()) {
            heap->insert( BranchSt(otherChild, new_distsq) );
        }

        /* Call recursively to search next level down. */
        searchLevel<with_removed>(result_set, vec, bestChild, mindist, checkCount, maxCheck, epsError, heap, checked);
    }

    /**
     * Performs an exact search in the tree starting from a node.
     */
    template<bool with_removed>
    void searchLevelExact(ResultSet<DistanceType>& result_set, const ElementType* vec, const NodePtr node, DistanceType mindist, const float epsError) const
    {
        /* If this is a leaf node, then do check and return. */
        if ((node->child1 == NULL)&&(node->child2 == NULL)) {
            int index = node->divfeat;
            if (with_removed) {
            	if (removed_points_.test(index)) return; // ignore removed points
            }
            DistanceType dist = distance_(node->point, vec, veclen_);
            result_set.addPoint(dist,index);

            return;
        }

        /* Which child branch should be taken first? */
        ElementType val = vec[node->divfeat];
        DistanceType diff = val - node->divval;
        NodePtr bestChild = (diff < 0) ? node->child1 : node->child2;
        NodePtr otherChild = (diff < 0) ? node->child2 : node->child1;

        /* Create a branch record for the branch not taken.  Add distance
            of this feature boundary (we don't attempt to correct for any
            use of this feature in a parent node, which is unlikely to
            happen and would have only a small effect).  Don't bother
            adding more branches to heap after halfway point, as cost of
            adding exceeds their value.
         */

        DistanceType new_distsq = mindist + distance_.accum_dist(val, node->divval, node->divfeat);

        /* Call recursively to search next level down. */
        searchLevelExact<with_removed>(result_set, vec, bestChild, mindist, epsError);

        if (mindist*epsError<=result_set.worstDist()) {
            searchLevelExact<with_removed>(result_set, vec, otherChild, new_distsq, epsError);
        }
    }
    
    void addPointToTree(NodePtr node, int ind)
    {
        ElementType* point = points_[ind];
        // std::cout << "pointer size: " << (points_[1]) - (points_[0]) << std::endl;
        // std::cout << "pointer 0: " << (points_[0]) << std::endl;
        // std::cout << "pointer 1: " << (points_[1]) << std::endl;
        
        if ((node->child1==NULL) && (node->child2==NULL)) {
            ElementType* leaf_point = node->point;
            ElementType max_span = 0;
            size_t div_feat = 0;
            for (size_t i=0;i<veclen_;++i) {
                ElementType span = std::abs(point[i]-leaf_point[i]);
                if (span > max_span) {
                    max_span = span;
                    div_feat = i;
                }
            }
            NodePtr left = new(pool_) Node();
            left->child1 = left->child2 = NULL;
            NodePtr right = new(pool_) Node();
            right->child1 = right->child2 = NULL;

            if (point[div_feat]<leaf_point[div_feat]) {
                left->divfeat = ind;
                left->point = point;
                left->point_num = ind;
                right->divfeat = node->divfeat;
                right->point = node->point;
                right->point_num = node->point_num;
            }
            else {
                left->divfeat = node->divfeat;
                left->point = node->point;
                left->point_num = node->point_num;
                right->divfeat = ind;
                right->point = point;
                right->point_num = ind;
            }
            node->divfeat = div_feat;
            node->divval = (point[div_feat]+leaf_point[div_feat])/2;
            node->child1 = left;
            node->child2 = right;            
        }
        else {
            if (point[node->divfeat]<node->divval) {
                addPointToTree(node->child1,ind);
            }
            else {
                addPointToTree(node->child2,ind);                
            }
        }
    }

    virtual void debug_index() {
        std::cout << "------------------------------" << std::endl;
        print_params(index_params_);
        std::cout << "------------------------------" << std::endl;

        std::cout << "number of flann datapoints: " << this->size_ << std::endl;
        std::cout << "size of one flann datapoint: " << this->veclen_ << std::endl;
        std::cout << "removed_count_: " << this->removed_count_ << std::endl;
        std::cout << "data_ptr_: " << this->data_ptr_ << std::endl;
        std::cout << "size of ids_: " << ids_.size() << std::endl;
        std::cout << "size of points_: " << points_.size() << std::endl;
        std::cout << "number of trees: " << trees_ << std::endl;
        std::cout << "number of tree roots: " << tree_roots_.size() << std::endl;
        
        // int tree_root_counter = 0;
        // for (auto it = begin (tree_roots_); it != end (tree_roots_); ++it) {
        //     std::cout << "tree root: " << tree_root_counter << std::endl;
        //     std::cout << std::addressof(*it) << std::endl; //print the address of tree root
        //     std::cout << (*it)->divfeat << std::endl; //print the dividing dimension of the node
		//     float point1 = (*it)->divval;
        //     std::cout << point1 << std::endl; //print the threshold value
                
        //     //In a kd tree, only the leaf nodes contain data.
        //     //The non-leaf nodes contain a dimension and threshold value (they do not contain data).
        //     int depth = 0;
        //     auto root = (*it);
        //     while (!(root->child1 == NULL && root->child2 == NULL)) { //traverse left until reaching a leaf node
        //         if (root->child1 != NULL) {
        //             root = root->child1;
        //             depth++;
        //         } else if (root->child2 != NULL) {
        //             root = root->child2;
        //             depth++;
        //         }
        //     }
        //     std::cout << "tree depth: " << depth << std::endl;
                
        //     tree_root_counter++;
        // }
    }

    /////////////////////////////////////////////////////////////
    //check_trees() - this function is used to validate that the
    //flattened tree matches the original tree
    void check_trees(int tree_root_num, char* ptr) {
        for(int i=0; i<2; i++) {
            //in-order tree traversal
            std::list<Node *> tree_nodes;
            Node *root;
            std::ofstream *outfile = new std::ofstream();

            if (i==0) {
                outfile->open("ref_tree.dat", std::ios::out | std::ios::binary | std::ios::trunc);
                root = tree_roots_[tree_root_num];
            } else {
                outfile->open("test_tree.dat", std::ios::out | std::ios::binary | std::ios::trunc);
                root = (Node*) ptr;
            }

            while (root != NULL || tree_nodes.size() != 0)
            {
                // Find the leftmost node
                while (root != NULL)
                {
                    tree_nodes.push_front(root); //inorder.push(root)
                    root = root->child1;         //root = root.left
                }

                root = tree_nodes.front();
                tree_nodes.pop_front(); //inorder.pop();

                if (root->child1 == NULL && root->child2 == NULL)
                {
                    //this is a leaf
                    outfile->write(reinterpret_cast<char *>(root->point), 256 * sizeof(float));
                }
                else
                {
                    //this is a non-leaf
                    outfile->write(reinterpret_cast<char *>(&root->divfeat), sizeof(int));
                    outfile->write(reinterpret_cast<char *>(&root->divval), sizeof(float));
                }

                root = root->child2; //root = root.right;
            }


            outfile->close();
        }
    }

    //#define STARTING_ADDR 0x6ffebb6f6000
       #define STARTING_ADDR 0x600000000000

    virtual void save_index(std::ofstream *outfile) 
    {
        //The save index file is in the following format:
        //---------------------------------------------------------------------------------------
        //| num visual words | visual words | num bytes tree 0 | pointer to tree 0 root | tree 0 bytes....
        //| (4 bytes)        | (many bytes) | (4 byte - int)   | (8 byte - pointer)     | (many bytes)
        //---------------------------------------------------------------------------------------
        //
        //---------------------------------------------------------------------------------------
        //| ....tree 0 bytes cont'd | num bytes tree 1 | pointer to tree 1 root | tree 1 bytes....
        //| (many bytes)            | (4 byte - int)   | (8 byte - int)         | (many bytes)
        //---------------------------------------------------------------------------------------
        //
        //There are 4 trees total and the root node addresses are stored in tree_roots_

        //Obtain the block of contiguous memory using mmap().  Malloc() does not return memory at a requested
        //address, but mmap() does.
        char *addr;  //pointer to the new memory
        unsigned long long int starting_addr = STARTING_ADDR;
        unsigned int length = 1300000000;  //1.3GB for now
        addr = (char *)mmap((void *)starting_addr, length, PROT_READ | PROT_WRITE, MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (addr == MAP_FAILED)
        {
            std::cout << "Error obtaining memory" << std::endl;
            exit(EXIT_FAILURE);
        } 
        else 
        {
            std::cout << "Successful memory mapping to: " << (unsigned long long) addr << std::endl;
        }

        //Copy the visual word data to the new memory.  Each visual word is a binary descriptor which consists of
        //256 floats, where each float is a binary number (0 or 1).  The visual words are copied into the new
        //memory first and then the trees are placed after the visual words.
        float *f_ptr, *f_write;
        f_write = (float*) addr;
        for(unsigned int i=0;i<points_.size();i++) 
        {
            f_ptr = (float*) points_[i];
            for (unsigned int j=0; j<256; j++) 
            {
                *f_write = f_ptr[j];
                f_write++;
            }
        }

        char *tree_addr_base = (char *)(addr + points_.size() * 256 * sizeof(float));
        char *tree_addr_cur, *region_start;
        tree_addr_cur = tree_addr_base;

        //map to store the mapping from the original memory address of a node, to the new memory address when it is loaded in again
        std::map<unsigned long long, unsigned long long> memory_addresses;
        struct Node *node_ptr;
        std::list<Node *> tree_nodes;  //used as a stack for in-order traversal of the tree
        Node *root;

        //This is the main loop that flattens the index trees.  The loop does an in-order traversal of 
        //the trees and copies the nodes to the memory block obtained using mmap().
        for (unsigned int tr = 0; tr<tree_roots_.size(); tr++)  //there are 4 trees
        {
            //at the memory region beginning, store the total number of bytes for the tree (4 byte int - does not include
            //the 8 byte pointer to the root node), then store the pointer to the root node of the tree (8 byte pointer)
            region_start = tree_addr_base; 
            tree_addr_base += sizeof(int) + sizeof(struct Node *); //make space to store the pointer to the root 
            tree_addr_cur = tree_addr_base;
            memory_addresses.clear();
            tree_nodes.clear();

            //in-order tree traversal
            root = tree_roots_[tr];
            int node_counter = 0;
            int leaf_counter = 0;

            while (root != NULL || tree_nodes.size() != 0)
            {
                // Find the leftmost node
                while (root != NULL)
                {
                    tree_nodes.push_front(root); //inorder.push(root)
                    root = root->child1;         //root = root.left
                }

                root = tree_nodes.front();
                tree_nodes.pop_front(); //inorder.pop();

                if (root->child1 == NULL && root->child2 == NULL)
                {
                    //this is a leaf
                    memory_addresses.insert(std::pair<unsigned long long, unsigned long long>((unsigned long long)root, (unsigned long long)tree_addr_cur));
                    node_ptr = (struct Node *)tree_addr_cur;
                    node_ptr->child1 = NULL;
                    node_ptr->child2 = NULL;
                    node_ptr->divfeat = root->divfeat;
                    node_ptr->divval = root->divval;
                    node_ptr->point_num = root->point_num;
                    node_ptr->point = (ElementType *)(starting_addr + (256 * sizeof(float) * node_ptr->point_num));
                    node_ptr->fix = 0;
                    leaf_counter++;
                    tree_addr_cur += sizeof(struct Node);
                }
                else
                {
                    //this is a non-leaf
                    // memory_addresses.insert(std::pair<unsigned long long,unsigned long long>((unsigned long long)root, (unsigned long long)tree_addr_cur));
                    memory_addresses[(unsigned long long)root] = (unsigned long long)tree_addr_cur;
                    node_ptr = (struct Node *)tree_addr_cur;

                    if (root->child1 != NULL)
                    {
                        //If this node has a left child, then that child has already been traversed.  Find
                        //the new address of the node and update the pointer.
                        node_ptr->child1 = (Node *)memory_addresses[(unsigned long long)root->child1];
                        memory_addresses.erase((unsigned long long)root->child1);
                    }
                    else
                    {
                        node_ptr->child1 = NULL;
                    }

                    if (root->child2 != NULL)
                    {
                        //If this node has a right child, then that child has not been traversed.  Set a flag 
                        //and placeholder address for updating in a final fixup pass through all the nodes.
                        node_ptr->child2 = root->child2; //have not traversed this node yet, so insert a placeholder to child2
                        node_ptr->fix = 1;               //mark this node for fixup later
                    }
                    else
                    {
                        node_ptr->child2 = NULL;
                        node_ptr->fix = 0;
                    }

                    node_ptr->divfeat = root->divfeat;
                    node_ptr->divval = root->divval;
                    node_ptr->point_num = root->point_num;
                    node_ptr->point = (ElementType *)NULL;
                    node_counter++;
                    tree_addr_cur += sizeof(struct Node);
                }

                root = root->child2; //root = root.right;
            }

            //Final loop pass to fix up the child2 pointers
            for (int i = 0; i < (node_counter + leaf_counter); i++)
            {
                Node *n_ptr = (Node *)(tree_addr_base + i * sizeof(struct Node));
                if (n_ptr->fix == 1)
                {
                    Node *new_ptr = (Node *)memory_addresses[(unsigned long long)n_ptr->child2];
                    memory_addresses.erase((unsigned long long)n_ptr->child2);
                    n_ptr->child2 = new_ptr;
                    n_ptr->fix = 0;
                }
            }

            //only the root node is left in the memory addresses map
            if (memory_addresses.size() != 1)
            {
                std::cout << "memory_addresses error: " << memory_addresses.size() << std::endl;
                exit(0);
            }

            int *size_ptr = (int*) region_start;
            *size_ptr = (node_counter + leaf_counter) * sizeof(struct Node); //store the number of bytes for the tree
            region_start += sizeof(int);
            unsigned long long *root_ptr = (unsigned long long *) region_start;
            std::map<unsigned long long, unsigned long long>::iterator it = memory_addresses.begin();
            *root_ptr = it->second;  //store the new address of the root node

            std::cout << "address of start of tree in memory: " << (unsigned long long)tree_addr_base << std::endl;
            std::cout << "address of root node: " << it->second << std::endl;
            std::cout << "address of end of tree in memory: " << (unsigned long long)tree_addr_cur << std::endl;
            if(tr ==0) {
                //check_trees(0, (char *)it->second);
            }

            std::cout << "nodes: " << node_counter << "  leaves: " << leaf_counter << std::endl;

            tree_addr_base = tree_addr_cur;
        }

        int total_flann_size = (unsigned long long)(tree_addr_cur - addr); //include visual words and 4 trees
        std::cout << "total flann size: " << total_flann_size << std::endl;

        //write out all the data to a file
        int point_size = points_.size();
        outfile->write(reinterpret_cast<char *>(&point_size), sizeof(int));
        outfile->write(reinterpret_cast<char *>(addr), total_flann_size);

        //write out the local variables of the class:  KDTreeIndex 
        outfile->write(reinterpret_cast<char *>(&trees_), sizeof(int));
        outfile->write(reinterpret_cast<char *>(&mean_), sizeof(DistanceType*));
        outfile->write(reinterpret_cast<char *>(&var_), sizeof(DistanceType*));

        //write out the local variables of the class:  NNIndex
        //Distance distance_;
        outfile->write(reinterpret_cast<char *>(&this->last_id_), sizeof(size_t));
        outfile->write(reinterpret_cast<char *>(&this->size_), sizeof(size_t));
        outfile->write(reinterpret_cast<char *>(&this->size_at_build_), sizeof(size_t));
        outfile->write(reinterpret_cast<char *>(&this->veclen_), sizeof(size_t));

        //IndexParams index_params_;
        // int indexparams_size = index_params_.size();
        // outfile->write(reinterpret_cast<char *>(&indexparams_size), sizeof(int));  //store the parameter string
        // for (std::map<std::string, any>::iterator iter = index_params_.begin(); iter != index_params_.end(); ++iter)
        // {
        //     int length = iter->first.length();
        //     outfile->write(reinterpret_cast<char *>(&length), sizeof(int));  //store the parameter string
        //     outfile->write(iter->first.data(), iter->first.length());  //store the parameter string
        //     outfile->write(reinterpret_cast<char *>(&iter->second), sizeof(any));  //store the value
        // }

        outfile->write(reinterpret_cast<char *>(&this->removed_), sizeof(bool));
        //DynamicBitset removed_points_;  //should be 0, so do not store
        outfile->write(reinterpret_cast<char *>(&this->removed_count_), sizeof(size_t));

        std::cout << "ids_ size: " << ids_.size() << std::endl;
        for (unsigned int i=0; i<ids_.size(); i++) {
            outfile->write(reinterpret_cast<char *>(&ids_[i]), sizeof(size_t));
        }

        outfile->write(reinterpret_cast<char *>(&this->data_ptr_), sizeof(ElementType*));

        //unmap the memory
        munmap(addr,length);

        // std::ofstream *idfile;
        // idfile = new std::ofstream();
        // idfile->open("ref_points.dat", std::ios::out | std::ios::binary | std::ios::trunc);
        // for (unsigned int i = 0; i < points_.size(); i++)
        // {
        //     idfile->write(reinterpret_cast<char *>(points_[i]), sizeof(float)*256);
        // }
        // idfile->close();
    }

    ///////////////////////////////////////////////////////////////////////
    //load_index
    virtual void load_index(std::ifstream *infile, char* data_ptr)
    {
        load_cached = 1;
        auto t1 = std::chrono::high_resolution_clock::now();

        //Call mmap() to have the memory block allocated at the same starting address.
        // char *addr;
        // unsigned long long int starting_addr = STARTING_ADDR;
        // unsigned int length = 1300000000;
        // addr = (char *) mmap((void *)starting_addr, length, PROT_READ | PROT_WRITE, MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        // if (addr == MAP_FAILED)
        // {
        //     std::cout << "Error" << std::endl;
        //     exit(EXIT_FAILURE);
        // } 
        // else 
        // {
        //     std::cout << "Successful mapping to: " << (unsigned long long) addr << std::endl;
        // }

        // //read in all the visual words from a file
        // int point_size;
        // infile->read(reinterpret_cast<char *>(&point_size), sizeof(int));

        // float* data_ptr = reinterpret_cast<float*> (addr);
        // for (int i=0; i<point_size; i++) 
        // {
        //     points_[i] = (ElementType*) data_ptr;
        //     infile->read(reinterpret_cast<char *>(data_ptr), 256 * sizeof(float));
        //     data_ptr += 256;
        // }


        // std::ofstream *idfile;
        // idfile = new std::ofstream();
        // idfile->open("test_points.dat", std::ios::out | std::ios::binary | std::ios::trunc);
        // for (int i = 0; i < point_size; i++)
        // {
        //     idfile->write(reinterpret_cast<char *>(points_[i]), sizeof(float)*256);
        // }
        // idfile->close();

        //read in the tree data for 4 trees
        std::cout << "--------size of trees_----------" << trees_ << std::endl;
        //this->tree_roots_.resize(trees_);
        this->tree_roots_.resize(trees_);
        NodePtr root[4] = {NULL, NULL, NULL, NULL};
        char* t_ptr = (char*) data_ptr;
        for (int i=0; i<4; i++) //number of trees
        { 
            int tree_size = 0;
            infile->read(reinterpret_cast<char *>(&tree_size), sizeof(int));  //read in the tree size in bytes
            std::cout << "loading tree size: " << tree_size << std::endl;
            infile->read(reinterpret_cast<char *>(&(root[i])), sizeof(NodePtr));  //read in the root node pointer address
            std::cout << "root address: " << root[i] << std::endl;
            this->tree_roots_.assign(i,root[i]);
            // tree_roots_[i] = root[i];

            t_ptr += sizeof(int) + sizeof(struct Node *);
            infile->read(reinterpret_cast<char *>(t_ptr), tree_size);  //read in the tree nodes in a contiguous block
            t_ptr += tree_size;
        }

        //write out the local variables of the class:  KDTreeIndex 
        infile->read(reinterpret_cast<char *>(&trees_), sizeof(int));
        infile->read(reinterpret_cast<char *>(&mean_), sizeof(DistanceType*));
        infile->read(reinterpret_cast<char *>(&var_), sizeof(DistanceType*));

        //write out the local variables of the class:  NNIndex
        //Distance distance_;
        infile->read(reinterpret_cast<char *>(&this->last_id_), sizeof(size_t));
        infile->read(reinterpret_cast<char *>(&this->size_), sizeof(size_t));
        infile->read(reinterpret_cast<char *>(&this->size_at_build_), sizeof(size_t));
        infile->read(reinterpret_cast<char *>(&this->veclen_), sizeof(size_t));

        //IndexParams index_params_;
        // int indexparams_size;
        // infile->read(reinterpret_cast<char *>(&indexparams_size), sizeof(int));  //restore the size of the parameters
        // for (int i=0; i<indexparams_size; i++)
        // {
        //     char buf[100];
        //     int length;
        //     infile->read(reinterpret_cast<char *>(&length), sizeof(int));  //read the length of the parameter string
        //     infile->read(reinterpret_cast<char *>(buf), length);  //read the parameter string
        //     buf[length] = 0; //NULL terminate the string;
        //     std::string s(buf);
        //     // any a;
        //     infile->read(reinterpret_cast<char *>(buf), sizeof(any));  //read the parameter value

        //     //TODO - index_params is already set on construction, so do not modify it
        //     // index_params_.insert(std::pair<std::string,any>(s,a));
        // }

        infile->read(reinterpret_cast<char *>(&this->removed_), sizeof(bool));
        //DynamicBitset removed_points_;  //should be 0, so do not store
        infile->read(reinterpret_cast<char *>(&this->removed_count_), sizeof(size_t));

        for (unsigned int i=0; i<ids_.size(); i++) {
            infile->read(reinterpret_cast<char *>(&ids_[i]), sizeof(size_t));
        }

        infile->read(reinterpret_cast<char *>(&this->data_ptr_), sizeof(ElementType*));

        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
        std::cout << "flann index load time in milliseconds: " << fp_ms.count() << std::endl;
    }

    virtual void set_cached (int cache) {
        this->cached = cache;
    }

private:
    void swap(KDTreeIndex& other)
    {
    	BaseClass::swap(other);
    	std::swap(trees_, other.trees_);
    	std::swap(tree_roots_, other.tree_roots_);
    	std::swap(pool_, other.pool_);
    }

private:

    enum
    {
        /**
         * To improve efficiency, only SAMPLE_MEAN random values are used to
         * compute the mean and variance at each level when building a tree.
         * A value of 100 seems to perform as well as using all values.
         */
        SAMPLE_MEAN = 100,
        /**
         * Top random dimensions to consider
         *
         * When creating random trees, the dimension on which to subdivide is
         * selected at random from among the top RAND_DIM dimensions with the
         * highest variance.  A value of 5 works well.
         */
        RAND_DIM=5
    };


    /**
     * Number of randomized trees that are used
     */
    int trees_;

    DistanceType* mean_;
    DistanceType* var_;

    /**
     * Array of k-d trees used to find neighbours.
     */
    std::vector<NodePtr> tree_roots_;

    /**
     * Pooled memory allocator.
     *
     * Using a pooled memory allocator is more efficient
     * than allocating memory directly when there is a large
     * number small of memory allocations.
     */
    PooledAllocator pool_;
    int load_cached;

    USING_BASECLASS_SYMBOLS
};   // class KDTreeIndex

}

#endif //FLANN_KDTREE_INDEX_H_
