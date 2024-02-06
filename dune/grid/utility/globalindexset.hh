// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 *
 *  \brief Provides a globally unique index for all entities of a distributed Dune grid
 *
 * Such functionality is relevant for a number of applications:
 *  - Map a degree of freedom associated with an entity to its
 *         location in a global matrix or global vector
 *  - Such indices for elements are needed as input to external mesh partitioners
 *  - Using matrix and vector routines from the PETSc or trilinos parallel linear algebra
 *    packages for distributed memory parallel computers.
 *
 *  Method: (1) The UniqueEntityPartition class assigns an owner process to each entity
 *
 *          (2) Compute the number of entities that are owned by each process
 *
 *          (3) we communicate the index of entities that are owned by the process to processes
 *              that also contain these entities but do not own them, so that on a non-owner process
 *              we have information on the index of the entity that it got from the owner-process;
 *
 *  \author    Benedikt Oswald, Patrick Leidenberger, Oliver Sander
 *
 *  \attention globally unique indices are ONLY provided for entities of the
 *             InteriorBorder_Partition type, NOT for the Ghost_Partition type !!!
 *
 *  \bug This code will most likely not work on grids with more than one element type.
 *       Please ask on the mailing list if you need this feature.
 *  \note The interface in this file is experimental, and may change without prior notice.
 */

#ifndef DUNE_GRID_UTILITY_GLOBALINDEXSET_HH
#define DUNE_GRID_UTILITY_GLOBALINDEXSET_HH

/** \brief Include standard header files. */
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <utility>
#include <algorithm>

/** include base class functionality for the communication interface */
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

/** include parallel capability */
#if HAVE_MPI
  #include <dune/common/parallel/mpihelper.hh>
#endif

namespace Dune
{

  /** \brief Calculate globally unique index over all processes in a Dune grid
   */
  template<class GridView>
  class GlobalIndexSet
  {
  public:
    /** \brief The number type used for global indices  */
    typedef int Index;

    /** \brief Helper class to provide access to subentity PartitionTypes with a run-time codimension
     *
     * This class can be removed if there is ever a method 'subPartitionType' similar to 'subIndex',
     * that takes a run-time codimension argument.
     */
    template <class Entity, int Codim>
    struct SubPartitionTypeProvider
    {
      /** \brief Get PartitionType of the i-th subentity of codimension 'codim' of entity 'entity'
       */
      static PartitionType get(const Entity& entity, int codim, int i)
      {
        if (codim==Codim)
          return entity.template subEntity<Codim>(i).partitionType();
        else
          return SubPartitionTypeProvider<Entity,Codim-1>::get(entity, codim, i);
      }
    };

    template <class Entity>
    struct SubPartitionTypeProvider<Entity,0>
    {
      static PartitionType get(const Entity& entity, int codim, int i)
      {
        return entity.template subEntity<0>(i).partitionType();
      }
    };

  private:
    /** define data types */
    typedef typename GridView::Grid Grid;

    typedef typename GridView::Grid::GlobalIdSet GlobalIdSet;
    typedef typename GridView::Grid::GlobalIdSet::IdType IdType;
    //typedef Dune::bigunsignedint<80> IdType;//std::size_t 
    typedef typename GridView::Traits::template Codim<0>::Iterator Iterator;

    typedef typename Grid::Communication Communication;

    typedef std::map<IdType,Index> MapId2Index;
    typedef std::map<Index,Index>    IndexMap;

    /*********************************************************************************************/
    /* calculate unique partitioning for all entities of a given codim in a given GridView,      */
    /* assuming they all have the same geometry, i.e. codim, type                                */
    /*********************************************************************************************/
    class UniqueEntityPartition
    {
    private:
      /* A DataHandle class to calculate the minimum of a std::vector which is accompanied by an index set */
      template<class IS, class V> // mapper type and vector type
      class MinimumExchange
      : public Dune::CommDataHandleIF<MinimumExchange<IS,V>,typename V::value_type>
      {
      public:
        //! export type of data for message buffer
        typedef typename V::value_type DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, unsigned int codim) const
        {
          return codim==indexSetCodim_;
        }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize (int dim, int codim) const
        {
          return true ;
        }

        /*! how many objects of type DataType have to be sent for a given entity
         *
         *  Note: Only the sender side needs to know this size. */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
          return 1 ;
        }

        /*! pack data from user to message buffer */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          buff.write(v_[indexset_.index(e)]);
        }

        /** \brief Unpack data from message buffer to user
         *
         * \param n The number of objects sent by the sender
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          DataType x;
          buff.read(x);
          if (x>=0) // other is -1 means, he does not want it
            v_[indexset_.index(e)] = std::min(x,v_[indexset_.index(e)]);
        }

        //! constructor
        MinimumExchange (const IS& indexset, V& v, unsigned int indexSetCodim)
        : indexset_(indexset),
          v_(v),
          indexSetCodim_(indexSetCodim)
        {}

      private:
        const IS& indexset_;
        V& v_;
        unsigned int indexSetCodim_;
      };

    public:
      /*! \brief Constructor needs to know the grid function space
       */
      UniqueEntityPartition (const GridView& gridview, unsigned int codim)
      : assignment_(gridview.size(codim))
      {
        /** extract types from the GridView data type */
        typedef typename GridView::IndexSet IndexSet;

        // assign own rank to entities that I might have
        for (auto it = gridview.template begin<0>(); it!=gridview.template end<0>(); ++it)
          for (unsigned int i=0; i<it->subEntities(codim); i++)
          {
            // Evil hack: I need to call subEntity, which needs the entity codimension as a static parameter.
            // However, we only have it as a run-time parameter.
            PartitionType subPartitionType = SubPartitionTypeProvider<typename GridView::template Codim<0>::Entity, GridView::dimension>::get(*it,codim,i);

            assignment_[gridview.indexSet().subIndex(*it,i,codim)]
              = ( subPartitionType==Dune::InteriorEntity or subPartitionType==Dune::BorderEntity )
              ? gridview.comm().rank()  // set to own rank
              : - 1;   // it is a ghost entity, I will not possibly own it.
          }

        /** exchange entity index through communication */
        MinimumExchange<IndexSet,std::vector<Index> > dh(gridview.indexSet(),assignment_,codim);

        gridview.communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);
      }

      /** \brief Which rank is the i-th entity assigned to? */
      int owner(size_t i)
      {
        return assignment_[i];
      }

      /** \brief Report the number of entities assigned to the rank 'rank' */
      size_t numOwners(int rank) const
      {
        return std::count(assignment_.begin(), assignment_.end(), rank);
      }

    private:
      std::vector<int> assignment_;
    };

  private:
    /* A DataHandle class to communicate the global index from the
     * owning to the non-owning entity; the class is based on the MinimumExchange
     * class in the parallelsolver.hh header file.
     */
    class IndexExchange
    : public Dune::CommDataHandleIF<IndexExchange,Index>
    {
    public:
      //! returns true if data for this codim should be communicated
      bool contains (int dim, unsigned int codim) const
      {
        return codim==indexSetCodim_;
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool fixedSize (int dim, int codim) const
      {
        return true;
      }

      /** \brief How many objects of type DataType have to be sent for a given entity
       *
       * \note Only the sender side needs to know this size.
       */
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        return 1;
      }

      /*! pack data from user to message buffer */
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        IdType id=globalidset_.id(e);

        if (indexSetCodim_==0)
          buff.write(mapid2entity_[id]);
        else
          buff.write((*mapid2entity_.find(id)).second);
      }

      /** \brief Unpack data from message buffer to user
       *
       * \param n The number of objects sent by the sender
       */
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& entity, size_t n)
      {
        Index x;
        buff.read(x);

        /** only if the incoming index is a valid one,
         *  i.e. if it is greater than zero, will it be
         *  inserted as the global index; it is made
         *  sure in the upper class, i.e. GlobalIndexSet,
         *  that non-owning processes use -1 to mark an entity
         *  that they do not own.
         */
        if(x >= 0) {
          const IdType id = globalidset_.id(entity);

          if (indexSetCodim_==0)
            mapid2entity_[id] = x;
          else
          {
            mapid2entity_.erase(id);
            mapid2entity_.insert(std::make_pair(id,x));

            const Index lindex = indexSet_.index(entity);
            localGlobalMap_[lindex] = x;
          }
        }
      }

      //! constructor
      IndexExchange (const GlobalIdSet& globalidset, MapId2Index& mapid2entity,
                     const typename GridView::IndexSet& localIndexSet, IndexMap& localGlobal,
                     unsigned int indexSetCodim)
      : globalidset_(globalidset),
      mapid2entity_(mapid2entity),
      indexSet_(localIndexSet),
      localGlobalMap_(localGlobal),
      indexSetCodim_(indexSetCodim)
      {}

    private:
      const GlobalIdSet& globalidset_;
      MapId2Index& mapid2entity_;

      const typename GridView::IndexSet& indexSet_;
      IndexMap& localGlobalMap_;
      unsigned int indexSetCodim_;
    };

  public:
    /** \brief Constructor for a given GridView
     *
     * This constructor calculates the complete set of global unique indices so that we can then
     *  later query the global index, by directly passing the entity in question.
     */
    GlobalIndexSet(const GridView& gridview, int codim, bool verbose = false)
    : gridview_(gridview),
      codim_(codim)
    {
          //bool verbose = false;
      int rank = gridview.comm().rank();
      int size = gridview.comm().size();

      //const typename GridView::IndexSet& indexSet = gridview.indexSet();

      std::unique_ptr<UniqueEntityPartition> uniqueEntityPartition;
      if (codim_!=0)
        uniqueEntityPartition = std::make_unique<UniqueEntityPartition>(gridview,codim_);//change from dune_release 2.7
        //uniqueEntityPartition = std::unique_ptr<UniqueEntityPartition>(new UniqueEntityPartition(gridview,codim_));

      int nLocalEntity = (codim_==0)
                    ? std::distance(gridview.template begin<0, Dune::Interior_Partition>(), gridview.template end<0, Dune::Interior_Partition>())
                    : uniqueEntityPartition->numOwners(rank);

      // Compute the global, non-redundant number of entities, i.e. the number of entities in the set
      // without double, aka. redundant entities, on the interprocessor boundary via global reduce. */
      nGlobalEntity_ = gridview.comm().template sum<int>(nLocalEntity);
	
        if (verbose)
          {
	  std::cout<<"check output start "<<rank<<" size "<<size<<" codim_ "<<codim_<<" nLocalEntity "<<nLocalEntity
	  <<" nGlobalEntity_ "<<nGlobalEntity_<<std::endl;
        }
		// for sharing the data
	std::vector<IdType>  sortedGlobalIdSet(nGlobalEntity_);//gather data on all threads
	std::vector<int>  P(size,0);//specifies the displacement (relative to recvbuf ) at which to place the incoming data from process i
      std::vector<int> offset(size);//to fill out P
      std::fill(offset.begin(), offset.end(), 0);

      /** Share number of locally owned entities */
      gridview_.comm().template allgather<int>(&nLocalEntity, 1, offset.data());

      for (int i=1; i<P.size(); i++)
	  {
		P[i] = P[i-1] +  offset[i-1];
	  }
	  int myoffset = P[rank];
          
        if (verbose)
          {
	  std::cout<<"checkP "<<rank<<" "<<myoffset<<" "<<P[rank]<<std::endl;
	  
      for (int i=0; i<P.size(); i++)
	  {
		  std::cout<<"P "<<P[i]<<" offset "<<offset[i]<<", ";
	  }
	  std::cout<<" offset[rank] "<<offset[rank]<<std::endl;
        }
		std::vector<IdType>  sendbuf(offset[rank]);//array with id of owned entities



      const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet();      /** retrieve globally unique Id set */

      Index localcontrib = 0;      /** initialize lid */

      if (codim_==0)  // This case is simpler
      {
        for (Iterator iter = gridview_.template begin<0>(); iter!=gridview_.template end<0>(); ++iter)
        {
          const IdType id = globalIdSet.id(*iter);      /** retrieve the entity's id */

          /** if the entity is owned by the process, go ahead with computing the global index */
          if (iter->partitionType() == Dune::InteriorEntity)
          {
        if (verbose)
          {
			  std::cout<<"sendbuf[localcontrib] = id "<<localcontrib<<" "<<id<<std::endl;
        }
			sendbuf[localcontrib] = id;
			localcontrib++;
          }

        }
      }
      else  // if (codim==0) else
      {
		  std::vector<bool> firstTime(gridview_.size(codim_));
		  std::fill(firstTime.begin(), firstTime.end(), true);

		  for(Iterator iter = gridview_.template begin<0>();iter!=gridview_.template end<0>(); ++iter)
		  {
        if (verbose)
          {
              std::cout<<"got iter"<<std::endl;
        }
			for (size_t i=0; i<iter->subEntities(codim_); i++)
			{
        if (verbose)
          {
              std::cout<<"got subEntities "<<i;
        }
			  IdType id=globalIdSet.subId(*iter,i,codim_);// real id

        if (verbose)
          {
                std::cout<<" id "<<id;
        }
			  Index idx = gridview_.indexSet().subIndex(*iter,i,codim_);/** local index?  */

        if (verbose)
          {
                std::cout<<", idx "<<idx;
        
                std::cout<<", sortedGlobalIdSet[idx] "<<sortedGlobalIdSet[idx];
        }
			  if (!firstTime[idx] )
				continue;

			  firstTime[idx] = false;

			  if (uniqueEntityPartition->owner(idx) == rank)  /** if the entity is owned by the process, go ahead */
			  {
				 const Index lindex = localcontrib;// idx;
                  
        if (verbose)
          {
                  std::cout<<" [sendbuf[lindex] = id "<<idx<<" "<<lindex<<" "<<id<<"] ";
        }
                    sendbuf[lindex] = id;
                    localcontrib++;
			  }
                
        if (verbose)
          {
            std::cout<<std::endl;
        }
			}

		  }
      }

        if (verbose)
          {
	  std::cout<<"checksendbuf "<<rank<<" "<<sendbuf.size()<<std::endl;
	  
      for (int i=0; i<sendbuf.size(); i++)
	  {
		  std::cout<<"sendbuf "<<sendbuf[i]<<", ";
	  }std::cout<<std::endl;
        }
	  
	  gridview_.comm().barrier();

      // 2nd stage of global index calculation: communicate global index for non-owned entities
		//int gatherv (const T* in, int sendlen, T* out, int* recvlen, int* displ) const
		gridview_.comm(). template allgatherv<IdType>(sendbuf.data(), nLocalEntity,sortedGlobalIdSet.data(),offset.data(), P.data());
          if (verbose)
          {
		std::cout<<"gatherv "<<rank<<" nLocalEntity "<<nLocalEntity<<" "<<P[rank]<<std::endl;
          }
		
		// Sort the vector in ascending order
		std::sort(sortedGlobalIdSet.begin(), sortedGlobalIdSet.end());
          if (verbose)
          {
	  std::cout<<"sortedGlobalIdSet "<<rank<<std::endl;
      for (int i=0; i<sortedGlobalIdSet.size(); i++)
	  {
		  std::cout<<"sortedGlobalIdSet "<<sortedGlobalIdSet[i]<<", ";
	  }std::cout<<"finished sorting the globalidset "<<std::endl;

          }
		// create and fill the map
      // initialize map that stores an entity's global index via it's globally unique id as key
      globalIndex_.clear();
		
      if (codim_==0)  // This case is simpler
      {
        	for(int contiguousGlobalIdSet = 0; contiguousGlobalIdSet <sortedGlobalIdSet.size(); contiguousGlobalIdSet++)
			{
				globalIndex_[sortedGlobalIdSet[contiguousGlobalIdSet]] = contiguousGlobalIdSet;
			}
      }
      else  // if (codim==0) else
      {
		  for(Iterator iter = gridview_.template begin<0>();iter!=gridview_.template end<0>(); ++iter)
		  {
              if (verbose)
          {
              std::cout<<"got iter"<<std::endl;
              }
			for (size_t i=0; i<iter->subEntities(codim_); i++)
			{
                if (verbose)
          {
              std::cout<<"got subEntities "<<i;
                }
			  IdType id=globalIdSet.subId(*iter,i,codim_);
                if (verbose)
          {
                std::cout<<" id "<<id;
                }
			  Index idx = gridview_.indexSet().subIndex(*iter,i,codim_);
                if (verbose)
          {
                std::cout<<", idx "<<idx;
                std::cout<<", sortedGlobalIdSet[idx] "<<sortedGlobalIdSet[idx];
                }

				globalIndex_.insert(std::make_pair(id,sortedGlobalIdSet[idx])); /** insert pair (key, value) into the map */

				const Index lindex = idx;
				localGlobalMap_[lindex] = sortedGlobalIdSet[idx];
                if (verbose)
          {
			  std::cout<<"globalIndex_["<<id<<"] "<<globalIndex_[id]<<", ";
			  std::cout<<"localGlobalMap_["<<idx<<"] "<<localGlobalMap_[idx]<<", ";

            std::cout<<std::endl;
                }
			}
		  }
      }
          if (verbose)
          {std::cout<<std::endl;
          std::cout<<"check outputend"<<std::endl;
          if (codim_==0)  // This case is simpler
          {
            for (Iterator iter = gridview_.template begin<0>(); iter!=gridview_.template end<0>(); ++iter)
            {
              const IdType id = globalIdSet.id(*iter);      /** retrieve the entity's id */
              std::cout<<"globalIndex_["<<id<<"] "<<globalIndex_[id]<<", ";
            }
          }
          else  // if (codim==0) else
          {
              for(Iterator iter = gridview_.template begin<0>();iter!=gridview_.template end<0>(); ++iter)
              {
                for (size_t i=0; i<iter->subEntities(codim_); i++)
                {
                  IdType id=globalIdSet.subId(*iter,i,codim_);

                  Index idx = gridview_.indexSet().subIndex(*iter,i,codim_);
                  std::cout<<"globalIndex_["<<id<<"] "<<globalIndex_[id]<<", ";
                  std::cout<<"localGlobalMap_["<<idx<<"] "<<localGlobalMap_[idx]<<", ";

                }

              }
          }std::cout<<std::endl;
          }
		
    }

    /** \brief Return the global index of a given entity */
    template <class Entity>
    Index index(const Entity& entity) const
    {
      if (codim_==0)
      {
        /** global unique index is only applicable for inter or border type entities */
        const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet(); /** retrieve globally unique Id set */
        const IdType id = globalIdSet.id(entity);                        /** obtain the entity's id */
        const Index gindex = globalIndex_.find(id)->second;                /** retrieve the global index in the map with the id as key */

        return gindex;
      }
      else
        return localGlobalMap_.find(gridview_.indexSet().index(entity))->second;
    }

    /** \brief Return the global index of a subentity of a given entity
     *
     * \param i Number of the requested subentity among all subentities of the given codimension
     * \param codim Codimension of the requested subentity
     */
    template <class Entity>
    Index subIndex(const Entity& entity, unsigned int i, unsigned int codim) const
    {
      if (codim_==0)
      {
        /** global unique index is only applicable for inter or border type entities */
        const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet(); /** retrieve globally unique Id set */
        const IdType id = globalIdSet.subId(entity,i,codim);                        /** obtain the entity's id */
        const Index gindex = globalIndex_.find(id)->second;                /** retrieve the global index in the map with the id as key */

        return gindex;
      }
      else
        return localGlobalMap_.find(gridview_.indexSet().subIndex(entity,i,codim))->second;
    }

    /** \brief Return the total number of entities over all processes that we have indices for
     *
     * \param codim If this matches GlobalIndexSet codimension, the number of entities is returned.
     *              Otherwise, zero is returned.
     */
    unsigned int size(unsigned int codim) const
    {
      return (codim_==codim) ? nGlobalEntity_ : 0;
    }
	
    // unsigned int getNglobalEntity() const
    // {
      // return nGlobalEntity_;
    // }
	
	// // to get localGlobalMap_[lindex] = sortedGlobalIdSet[idx];
	// IndexMap getLocalGlobalMap() const
	// {
		// return localGlobalMap_;
	// }

	// // globalIndex_.insert(std::make_pair(id,sortedGlobalIdSet[idx]));
	// //MapId2Index getGlobalIndex() const
	// MapId2Index getGlobalIndex() const 
	// {
		// return globalIndex_;
	// }
	
  protected:
    const GridView gridview_;

    /** \brief Codimension of the entities that we hold indices for */
    unsigned int codim_;

    //! Global number of entities, i.e. number of entities without rendundant entities on interprocessor boundaries
    int nGlobalEntity_;

    IndexMap localGlobalMap_;

    /** \brief Stores global index of entities with entity's globally unique id as key
     */
    //std::map<Dune::bigunsignedint<80>,int>   globalIndex_;
	MapId2Index globalIndex_;
  };

}  // namespace Dune

#endif /* DUNE_GRID_UTILITY_GLOBALINDEXSET_HH */
