#ifndef DUNE_SPGRID_LINKAGE_HH
#define DUNE_SPGRID_LINKAGE_HH

#include <vector>

#include <dune/grid/spgrid/partitionlist.hh>
#include <dune/grid/spgrid/partitionpool.hh>

namespace Dune
{

  // SPLinkage
  // ---------

  template< int dim >
  class SPLinkage
  {
    typedef SPLinkage< dim > This;

  public:
    typedef SPPartitionPool< dim > PartitionPool;
    typedef SPPartitionList< dim > PartitionList;

    typedef typename PartitionPool::Mesh Mesh;
    typedef typename PartitionPool::MultiIndex MultiIndex;

    class Interface;

    SPLinkage ( const int localRank,
                const PartitionPool &localPool,
                const std::vector< Mesh > &decomposition );

    const Interface &interface ( const InterfaceType iftype ) const;

  private:
    template< InterfaceType iftype >
    bool build ( const int localRank, const PartitionPool &localPool,
                 const int remoteRank, const PartitionPool &remotePool );

   const PartitionList *
   intersect ( const bool order, const PartitionList &local, const PartitionList &remote ) const;

    // note: We use the knowledge that interfaces are numbered 0, ..., 4.
    Interface interface_[ 5 ];
  };



  // SPLinkage::Interface
  // --------------------

  template< int dim >
  class SPLinkage< dim >::Interface
  {
    struct Node;

    typedef std::vector< Node > NodeContainer;

  public:
    typedef typename NodeContainer::const_iterator Iterator;

    Interface ();
    Interface ( const Interface &other );

    ~Interface ();

    Iterator begin () const { return nodes_.begin(); }
    Iterator end () const { return nodes_.end(); }

    std::size_t size () const { return nodes_.size(); }

    void add ( int rank, const PartitionList *sendList, const PartitionList *receiveList )
    {
      nodes_.emplace_back( rank, sendList, receiveList );
    }

  private:
    NodeContainer nodes_;
  };



  // SPLinkage::Interface::Node
  // --------------------------

  template< int dim >
  struct SPLinkage< dim >::Interface::Node
  {
    static_assert( (int( ForwardCommunication ) == 0) && (int( BackwardCommunication ) == 1),
                   "enumeration CommunicationDirection has changed." );

    Node ( const int rank, const PartitionList *sendList, const PartitionList *receiveList );

    int rank () const;
    const PartitionList &sendList ( const CommunicationDirection direction = ForwardCommunication ) const;
    const PartitionList &receiveList ( const CommunicationDirection direction = ForwardCommunication ) const;

    void destroy ();

  private:
    int rank_;
    const PartitionList *partitionList_[ 2 ];
  };



  // SPCommunicationInterface
  // ------------------------

  template< InterfaceType iftype >
  struct SPCommunicationInterface;

  template<>
  struct SPCommunicationInterface< InteriorBorder_InteriorBorder_Interface >
  {
    static const PartitionIteratorType sendPartition = InteriorBorder_Partition;
    static const PartitionIteratorType receivePartition = InteriorBorder_Partition;
  };

  template<>
  struct SPCommunicationInterface< InteriorBorder_All_Interface >
  {
    static const PartitionIteratorType sendPartition = InteriorBorder_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };

  template<>
  struct SPCommunicationInterface< Overlap_OverlapFront_Interface >
  {
    static const PartitionIteratorType sendPartition = Overlap_Partition;
    static const PartitionIteratorType receivePartition = OverlapFront_Partition;
  };

  template<>
  struct SPCommunicationInterface< Overlap_All_Interface >
  {
    static const PartitionIteratorType sendPartition = Overlap_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };

  template<>
  struct SPCommunicationInterface< All_All_Interface >
  {
    static const PartitionIteratorType sendPartition = All_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };



  // Implementation of SPLinkage
  // ---------------------------

  template< int dim >
  inline SPLinkage< dim >
    ::SPLinkage ( const int localRank,
                  const PartitionPool &localPool,
                  const std::vector< Mesh > &decomposition )
  {
    const int size = decomposition.size();

    const Mesh &globalMesh = localPool.globalMesh();
    const MultiIndex &overlap = localPool.overlap();

    for( int remoteRank = 0; remoteRank < size; ++remoteRank )
    {
      if( remoteRank == localRank )
        continue;
      PartitionPool remotePool( decomposition[ remoteRank ], globalMesh, overlap, localPool.topology() );
      if( build< All_All_Interface >( localRank, localPool, remoteRank, remotePool ) )
      {
        build< InteriorBorder_InteriorBorder_Interface >( localRank, localPool, remoteRank, remotePool );
        build< InteriorBorder_All_Interface >( localRank, localPool, remoteRank, remotePool );
        build< Overlap_OverlapFront_Interface >( localRank, localPool, remoteRank, remotePool );
        build< Overlap_All_Interface >( localRank, localPool, remoteRank, remotePool );
      }
    }
  }


  template< int dim >
  inline const typename SPLinkage< dim >::Interface &
  SPLinkage< dim >::interface ( const InterfaceType iftype ) const
  {
    assert( (int( iftype ) >= 0) && (int( iftype ) < 5) );
    return interface_[ int( iftype ) ];
  }


  template< int dim >
  template< InterfaceType iftype >
  inline bool SPLinkage< dim >
    ::build ( const int localRank, const PartitionPool &localPool,
              const int remoteRank, const PartitionPool &remotePool )
  {
    const PartitionIteratorType piSend = SPCommunicationInterface< iftype >::sendPartition;
    const PartitionIteratorType piReceive = SPCommunicationInterface< iftype >::receivePartition;

    const bool order = (localRank < remoteRank);

    // build intersection lists
    const PartitionList *sendList, *receiveList;
    sendList = intersect( order, localPool.template get< piSend >(), remotePool.template get< piReceive >() );
    if( piSend != piReceive )
      receiveList = intersect( order, localPool.template get< piReceive >(), remotePool.template get< piSend >() );
    else
      receiveList = sendList;

    // if both lists are empty, no communication is necessary
    if( sendList->empty() && receiveList->empty() )
    {
      if( sendList != receiveList )
        delete receiveList;
      delete sendList;
      return false;
    }

    assert( (iftype >= 0) && (iftype < 5) );
    interface_[ iftype ].add( remoteRank, sendList, receiveList );
    return true;
  }


  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList *
  SPLinkage< dim >::intersect ( const bool order, const PartitionList &local, const PartitionList &remote ) const
  {
    typedef SPBasicPartition< dim > Intersection;
    typedef typename PartitionList::Iterator Iterator;
    typedef typename PartitionList::Partition Partition;

    // order = true <=>  pit local iterator, qit remote iterator

    PartitionList *link = new PartitionList;
    for( Iterator pit = (order ? local.begin() : remote.begin()); pit; ++pit )
    {
      for( Iterator qit = (order ? remote.begin() : local.begin()); qit; ++qit )
      {
        const int number = (order ? pit->number() : qit->number());
        Intersection intersection = pit->intersect( *qit );
        if( !intersection.empty() )
          *link += Partition( intersection, number );
      }
    }
    return link;
  }



  // Implementation of SPLinkage::Interface
  // --------------------------------------

  template< int dim >
  inline SPLinkage< dim >::Interface::Interface ()
  {}


  template< int dim >
  inline SPLinkage< dim >::Interface::Interface ( const Interface &other )
  {
    nodes_.reserve( other.nodes_.size() );
    const Iterator end = other.end();
    for( Iterator it = other.begin(); it != end; ++it )
    {
      const PartitionList *sendList = new PartitionList( it->sendList() );
      const PartitionList *receiveList = new PartitionList( it->receiveList() );
      nodes_->push_back( it->rank(), sendList, receiveList );
    }
  }


  template< int dim >
  inline SPLinkage< dim >::Interface::~Interface ()
  {
    typedef typename NodeContainer::iterator Iterator;
    const Iterator end = nodes_.end();
    for( Iterator it = nodes_.begin(); it != end; ++it )
      it->destroy();
  }



  // Implementation of SPLinkage::Interface
  // --------------------------------------

  template< int dim >
  inline SPLinkage< dim >::Interface::Node
    ::Node ( const int rank, const PartitionList *sendList, const PartitionList *receiveList )
  : rank_( rank )
  {
    partitionList_[ 0 ] = sendList;
    partitionList_[ 1 ] = receiveList;
  }


  template< int dim >
  inline int SPLinkage< dim >::Interface::Node::rank () const
  {
    assert( rank_ >= 0 );
    return rank_;
  }


  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList &
  SPLinkage< dim >::Interface::Node::sendList ( const CommunicationDirection direction ) const
  {
    assert( (int( direction ) >= 0) && (int( direction ) <= 1) );
    assert( partitionList_[ direction ] != 0 );
    return *partitionList_[ direction ];
  }

  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList &
  SPLinkage< dim >::Interface::Node::receiveList ( const CommunicationDirection direction ) const
  {
    assert( (int( direction ) >= 0) && (int( direction ) <= 1) );
    assert( partitionList_[ 1 - direction ] != 0 );
    return *partitionList_[ 1 - direction ];
  }


  template< int dim >
  inline void SPLinkage< dim >::Interface::Node::destroy ()
  {
    rank_ = -1;
    if( partitionList_[ 0 ] != partitionList_[ 1 ] )
      delete partitionList_[ 1 ];
    delete partitionList_[ 0 ];
    partitionList_[ 0 ] = partitionList_[ 1 ] = 0;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_LINKAGE_HH
