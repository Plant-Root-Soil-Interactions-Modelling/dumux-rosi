#ifndef DUNE_SPGRID_PARTITIONLIST_HH
#define DUNE_SPGRID_PARTITIONLIST_HH

#include <iterator>

#include <dune/grid/spgrid/partition.hh>

namespace Dune
{

  // SPPartitionList
  // ---------------

  template< int dim >
  class SPPartitionList
  {
    typedef SPPartitionList< dim > This;

  protected:
    struct Node;

  public:
    typedef SPPartition< dim > Partition;

    typedef typename Partition::MultiIndex MultiIndex;
    typedef typename Partition::Mesh Mesh;
    
    struct Iterator;

    SPPartitionList () : head_( nullptr ) {}

    SPPartitionList ( const This &other ) : head_( other.head_ ? new Node( *other.head_ ) : nullptr ) {}

    SPPartitionList ( This &&other ) : head_( other.head_ ) { other.head_ = nullptr; }

    ~SPPartitionList () { delete head_; }

    This &operator= ( const This &other )
    {
      delete head_;
      head_ = (other.head_ ? new Node( *other.head_ ) : nullptr);
      return *this;
    }

    This &operator= ( This &&other )
    {
      delete head_;
      head_ = other.head_;
      other.head_ = nullptr;
      return *this;
    }

    This &operator+= ( const Partition &partition );

    Iterator begin () const { return Iterator( head_ ); }
    Iterator end () const { return Iterator( nullptr ); }

    bool contains ( const MultiIndex &id, unsigned int number ) const;
    const Partition *findPartition ( const MultiIndex &id ) const;
    int volume () const;

    bool empty () const { return !head_; }
    unsigned int size () const;

  protected:
    Node *head_;
  };



  // SPPartitionList::Node
  // ---------------------

  template< int dim >
  struct SPPartitionList< dim >::Node
  {
    explicit Node ( const Partition &partition )
      : partition_( partition ),
        next_( nullptr )
    {}

    Node ( const Node &other )
      : partition_( other.partition_ ),
        next_( other.next_ ? new Node( *other.next_ ) : nullptr )
    {}

    Node ( Node &&other )
      : partition_( other.partition_ ),
        next_( other.next_ )
    {
      other.next_ = nullptr;
    }

    ~Node () { delete next_; }

    Node &operator= ( const Node & ) = delete;
    Node &operator= ( Node && ) = delete;

    void append ( Node *other )
    {
      if( next_ )
        next_->append( other );
      else
        next_ = other;
    }

    const Partition &partition () const { return partition_; }

    const Node *next () const { return next_; }

  private:
    Partition partition_;
    Node *next_;
  };



  // SPPartitionList::Iterator
  // -------------------------

  template< int dim >
  struct SPPartitionList< dim >::Iterator
    //: public std::iterator< std::forward_iterator_tag, const Partition >
  {
	  
    using iterator_category = std::forward_iterator_tag;
    using value_type = const Partition;
    using difference_type = std::ptrdiff_t;
    using pointer = const Partition *;
    using reference = const Partition &;

    explicit Iterator ( const Node *node = nullptr )
      : node_( node )
    {}

    Iterator &operator++ ()
    {
      assert( *this );
      node_ = node_->next();
      return *this;
    }

    Iterator operator++ ( int ) { Iterator copy( *this ); ++(*this); return copy; }

    operator bool () const { return bool( node_ ); }

    bool operator== ( const Iterator &other ) const { return (node_ == other.node_); }
    bool operator!= ( const Iterator &other ) const { return (node_ != other.node_); }

    const Partition &operator* () const { assert( *this ); return node_->partition(); }
    const Partition *operator-> () const { assert( *this ); return &(node_->partition()); }

  private:
    const Node *node_;
  };



  // Implementation of SPPartitionList
  // ---------------------------------

  template< int dim >
  inline typename SPPartitionList< dim >::This &
  SPPartitionList< dim >::operator+= ( const Partition &partition )
  {
    if( head_ )
      head_->append( new Node( partition ) );
    else
      head_ = new Node( partition );
    return *this;
  }


  template< int dim >
  inline bool
  SPPartitionList< dim >
    ::contains ( const MultiIndex &id, unsigned int number ) const
  {
    const Partition *partition = findPartition( id );
    assert( !partition || (partition->number() == number) );
    return bool( partition );
  }


  template< int dim >
  inline const typename SPPartitionList< dim >::Partition *
  SPPartitionList< dim >::findPartition ( const MultiIndex &id ) const
  {
    for( const Node *it = head_; it; it = it->next() )
    {
      if( it->partition().contains( id ) )
        return &(it->partition());
    }
    return nullptr;
  }


  template< int dim >
  inline int SPPartitionList< dim >::volume () const
  {
    int volume = 0;
    for( const Node *it = head_; it; it = it->next() )
      volume += it->partition().volume();
    return volume;
  }


  template< int dim >
  inline unsigned int SPPartitionList< dim >::size () const
  {
    unsigned int size = 0;
    for( const Node *it = head_; it; it = it->next() )
      ++size;
    return size;
  }



  // Auxilliary Functions for SPPartitionList
  // ----------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPPartitionList< dim > &partition )
  {
    typedef typename SPPartitionList< dim >::Iterator Iterator;
    std::string separator = "";
    for( Iterator it = partition.begin(); it; ++it )
    {
      out << separator << *it;
      separator = "; ";
    }
    return out;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
