#ifndef DUMUX_MPFA_O_INTERACTIONVOLUMECONTAINER_HH
#define DUMUX_MPFA_O_INTERACTIONVOLUMECONTAINER_HH


#include <dumux/common/math.hh>
#include <dumux/implicit/mpfa/mpfaproperties.hh>

namespace Dumux{

template<class TypeTag>
class MpfaOInteractionVolumeContainer
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolume) BoundaryInteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolumeManager) InteractionVolumeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolumeManager) BoundaryInteractionVolumeManager;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryLayers) BoundaryLayers;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum
    {
        dim = GridView::dimension,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

public:

    MpfaOInteractionVolumeContainer(): problemPtr_(0)
    {}

    void init(Problem &problem)
    {
        problemPtr_ = &problem;
        innerInteractionVolumes_.clear();
        boundaryInteractionVolumes_.clear();
        innerVolumesMap_.clear();
        boundaryVolumesMap_.clear();

        int innerVolumesCounter = 0;
        int boundaryVolumesCounter = 0;

        // first of all, do only the elements, that have intersections with the boundary
        // in order to properly set up the interaction volumes on the boundary
        ElementIterator eIt = problem.gridView().template begin<0>();
        ElementIterator eEndIt = problem.gridView().template end<0>();
        for(; eIt != eEndIt; ++eIt)
        {
            if (!eIt->hasBoundaryIntersections())
                continue;
            // Reference element will be needed later on
            const ReferenceElement& referenceElement = ReferenceElements::general(eIt->geometry().type());

            IntersectionIterator isIt = problem.gridView().ibegin(*eIt);
            IntersectionIterator isEnd = problem.gridView().iend(*eIt);
            for(; isIt != isEnd; ++isIt)
            {
                // loop over nodes of facet
                for(int nodeIndex = 0; nodeIndex < isIt->geometry().corners(); nodeIndex++)
                {
                    // get node index w.r.t element
                    int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, nodeIndex, dim);

                    // get global index of central vertex
                    int centralVertGlobalIdx = problem_().vertexMapper().index( eIt->template subEntity<dim>(localVertCorner) );

                    if ( hasVolumeBeenStored(centralVertGlobalIdx) )
                        continue;

                    // If node is on the boundary, make boundary interaction volumes for
                    // the different phases
                    if(nodeTouchesBoundary(*eIt, isIt, localVertCorner))
                    {
                        Dune::FieldVector<std::shared_ptr<BoundaryInteractionVolume>, numPhases> boundaryVolumes;
                        for(int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
                        {
                            boundaryVolumes[phaseIdx] = std::make_shared<BoundaryInteractionVolume>();
                            BoundaryInteractionVolumeManager::fillInteractionVolume(boundaryVolumes[phaseIdx], *eIt, problem, isIt->indexInInside(), localVertCorner, centralVertGlobalIdx, phaseIdx);
                        }
                        boundaryInteractionVolumes_.push_back(boundaryVolumes);
                        boundaryVolumesMap_.insert(std::pair<int, int> (centralVertGlobalIdx, boundaryVolumesCounter));
                        boundaryVolumesCounter++;
                    }
                    else
                    {
                        std::shared_ptr<InteractionVolume> interActionVolume = std::make_shared<InteractionVolume>();
                        InteractionVolumeManager::fillInteractionVolume(interActionVolume, *eIt, problem, isIt->indexInInside(), localVertCorner, centralVertGlobalIdx, /*dummy*/0);
                        innerInteractionVolumes_.push_back(interActionVolume);
                        innerVolumesMap_.insert(std::pair<int, int> (centralVertGlobalIdx, innerVolumesCounter));
                        innerVolumesCounter++;
                    }
                }
            }
        }

        // Now set up the remaining inner interaction volumes
        eIt = problem.gridView().template begin<0>();
        eEndIt = problem.gridView().template end<0>();
        for(; eIt != eEndIt; ++eIt)
        {
            if (eIt->hasBoundaryIntersections())
                continue;

            // Reference element will be needed later on
            const ReferenceElement& referenceElement = ReferenceElements::general(eIt->geometry().type());

            IntersectionIterator isIt = problem.gridView().ibegin(*eIt);
            IntersectionIterator isEnd = problem.gridView().iend(*eIt);
            for(; isIt != isEnd; ++isIt)
            {
                // loop over nodes of facet
                for(int nodeIndex = 0; nodeIndex < isIt->geometry().corners(); nodeIndex++)
                {
                    // get node index w.r.t element
                    int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, nodeIndex, dim);

                    // get global index of central vertex
                    int centralVertGlobalIdx = problem_().vertexMapper().index( eIt->template subEntity<dim>(localVertCorner) );

                    if ( hasVolumeBeenStored(centralVertGlobalIdx) )
                        continue;

                    std::shared_ptr<InteractionVolume> interActionVolume = std::make_shared<InteractionVolume>();
                    InteractionVolumeManager::fillInteractionVolume(interActionVolume, *eIt, problem, isIt->indexInInside(), localVertCorner, centralVertGlobalIdx, /*dummy*/0);
                    innerInteractionVolumes_.push_back(interActionVolume);
                    innerVolumesMap_.insert(std::pair<int, int> (centralVertGlobalIdx, innerVolumesCounter));
                    innerVolumesCounter++;
                }
            }
        }
    }

    bool nodeTouchesBoundary(const Element &element, IntersectionIterator &isIt, int localNodeIdx)
    {
        if (isIt->boundary())
            return true;
        else
        {
            // the reference element
            const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

            IntersectionIterator isIt2 = problem_().gridView().ibegin(element);
            IntersectionIterator isIt2End = problem_().gridView().iend(element);
            for(; isIt2 != isIt2End; ++isIt2)
            {
                if (isIt2 == isIt)
                    continue;
                else
                    for(int corner = 0; corner < isIt2->geometry().corners(); corner++)
                    {
                        // get node index w.r.t element
                        int localVertCorner = referenceElement.subEntity(isIt2->indexInInside(), 1, corner, dim);

                        if (localVertCorner == localNodeIdx && !isIt2->neighbor())
                            return true;
                    }
            }
        }

        return false;
    }

    int isBoundaryVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        // check the boundary volumes for the globalVertexIndex
        std::map<int, int>::const_iterator boundaryIterator = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIterator != boundaryVolumesMap_.end())
            return BoundaryLayers::boundary;
        else
            return BoundaryLayers::interior;
    }

    int isBoundaryVolume(int globalVertIdx) const
    {
        // check the boundary volumes for the globalVertexIndex
        std::map<int, int>::const_iterator boundaryIterator = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIterator != boundaryVolumesMap_.end())
            return BoundaryLayers::boundary;
        else
            return BoundaryLayers::interior;
    }

    // Gives back an interaction volume which is inside the domain
    std::shared_ptr<InteractionVolume> getInsideInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );
        std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
        if(innerIt != innerVolumesMap_.end())
            return innerInteractionVolumes_[innerIt->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find inner interaction volume!!!");
    }
    // This function is identical to the one above, but is needed for compatibility with other mpfa methods
    std::shared_ptr<InteractionVolume> getInnerInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );
        std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
        if(innerIt != innerVolumesMap_.end())
            return innerInteractionVolumes_[innerIt->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find inner interaction volume!!!");
    }
    std::shared_ptr<InteractionVolume> getInsideInteractionVolumeFromGlobalVertexIndex(int globalVertIdx) const
    {
        std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
        if(innerIt != innerVolumesMap_.end())
            return innerInteractionVolumes_[innerIt->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find inner interaction volume!!!");
    }

    // Gives back a boundary interaction volume for the respective equation
    std::shared_ptr<BoundaryInteractionVolume> getBoundaryInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx, int phaseIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIt != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[boundaryIt->second][phaseIdx];
        else
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find boundary interaction volume!!!");
    }
    std::shared_ptr<BoundaryInteractionVolume> getBoundaryInteractionVolumeFromGlobalVertexIndex(int globalVertIdx, int phaseIdx) const
    {
        std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIt != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[boundaryIt->second][phaseIdx];
        else
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find boundary interaction volume from vertex!!!");
    }

    // Gives back any interaction volume in the domain, independent of its location and of the phase- or eqIdx
    // This method is necessary for getting informations on the stencil anywhere in the domain
    std::shared_ptr<BoundaryInteractionVolume> getMpfaOInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIt != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[boundaryIt->second][0];
        else
        {
            std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
            if(innerIt != innerVolumesMap_.end())
                return innerInteractionVolumes_[innerIt->second];
            else
                DUNE_THROW(Dune::InvalidStateException, "Couldn't find inner interaction volume!!!");
        }
    }
    std::shared_ptr<BoundaryInteractionVolume> getMpfaOInteractionVolumeFromGlobalVertexIndex(int globalVertIdx) const
    {
        std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
        if(boundaryIt != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[boundaryIt->second][0];
        else
        {
            std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
            if(innerIt != innerVolumesMap_.end())
                return innerInteractionVolumes_[innerIt->second];
            else
                DUNE_THROW(Dune::InvalidStateException, "Couldn't find inner interaction volume!!!");
        }
    }

    bool hasVolumeBeenStored(int globalVertIdx) const
    {
        // first check the inner vertices, then the boundary vertices
        std::map<int, int>::const_iterator innerIt = innerVolumesMap_.find(globalVertIdx);
        if (innerIt != innerVolumesMap_.end())
            return true;
        else
        {
            std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
            if(boundaryIt != boundaryVolumesMap_.end())
                return true;
            else
                return false;
        }
    }

protected:
    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

private:
    Problem *problemPtr_;
    std::vector< std::shared_ptr<InteractionVolume> > innerInteractionVolumes_;
    std::vector< Dune::FieldVector<std::shared_ptr<BoundaryInteractionVolume>, numPhases> > boundaryInteractionVolumes_;

    std::map<int, int> innerVolumesMap_;
    std::map<int, int> boundaryVolumesMap_;
};
}
#endif /* DUMUX_MPFA_O_INTERACTIONVOLUMECONTAINER_HH */
