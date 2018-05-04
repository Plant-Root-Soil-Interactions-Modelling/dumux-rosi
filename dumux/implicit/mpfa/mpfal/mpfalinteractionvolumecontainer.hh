#ifndef DUMUX_MPFA_L_INTERACTIONVOLUMECONTAINER_HH
#define DUMUX_MPFA_L_INTERACTIONVOLUMECONTAINER_HH


#include <dumux/common/math.hh>
#include <dumux/implicit/mpfa/mpfaproperties.hh>

#include <dumux/implicit/mpfa/mpfao/2d/mpfao2dinteractionvolume.hh>
#include <dumux/implicit/mpfa/mpfao/2d/mpfao2dmanager.hh>

namespace Dumux{

template<class TypeTag>
class MpfaLInteractionVolumeContainer
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolume) BoundaryInteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolumeManager) InteractionVolumeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolumeManager) BoundaryManager;
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

    MpfaLInteractionVolumeContainer(): problemPtr_(0)
    {}

    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        // For the Mpfa-L method we use O-Method interaction volumes for the sub faces
        // that are connected to the boundary and L-Method interaction volumes for interior faces
        int boundaryVolumesCounter = 0;
        int innerVolumesCounter = 0;
        int intermediateVolumesCounter = 0;
        innerInteractionVolumes_.clear();
        boundaryInteractionVolumes_.clear();
        intermediateInteractionVolumes_.clear();
        innerVolumesMap_.clear();
        boundaryVolumesMap_.clear();
        intermediateVolumesMap_.clear();

        // first of all, do only the elements, that have intersections with the boundary
        // in order to properly set up MPFA-O interaction volumes
        ElementIterator eIt = problem.gridView().template begin<0>();
        ElementIterator eEndIt = problem.gridView().template end<0>();
        for(; eIt != eEndIt; ++eIt)
        {
            if (!eIt->hasBoundaryIntersections())
                continue;

            // global element index
            int globalElemIdx = problem_().elementMapper().index(*eIt);

            // Reference element will be needed later on
            const ReferenceElement& referenceElement = ReferenceElements::general(eIt->geometry().type());

            IntersectionIterator isIt = problem.gridView().ibegin(*eIt);
            IntersectionIterator isEnd = problem.gridView().iend(*eIt);
            for(; isIt != isEnd; ++isIt)
            {
                // loop over the nodes of the facet
                for(int facetNode = 0; facetNode < isIt->geometry().corners(); facetNode++)
                {
                    // get node index w.r.t element
                    int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, facetNode, dim);

                    // get global index
                    int globalVertIdx = problem_().vertexMapper().index( eIt->template subEntity<dim>(localVertCorner) );

                    // only make volume if does not already exist
                    if (!hasVolumeBeenStored(*eIt, isIt->indexInInside(), localVertCorner) )
                    {
                        // If node touches boundary, create a vector of boundary interaction volume
                        if ( nodeTouchesBoundary(*eIt, isIt, localVertCorner) )
                        {
                            Dune::FieldVector< std::shared_ptr<BoundaryInteractionVolume>, numPhases> boundaryVolumes;
                            for(int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
                            {
                                boundaryVolumes[phaseIdx] = std::make_shared<BoundaryInteractionVolume>();
                                BoundaryManager::fillInteractionVolume(boundaryVolumes[phaseIdx], *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx, phaseIdx);
                            }
                            boundaryInteractionVolumes_.push_back(boundaryVolumes);
                            boundaryVolumesMap_.insert( std::pair<int, int> (globalVertIdx, boundaryVolumesCounter) );
                            boundaryVolumesCounter++;
                        }
                        // if the node belongs to influence region of the boundary, create boundary interaction volume
                        else if ( GET_PROP_VALUE(TypeTag, FullMpfaOBoundary) && fullFaceTouchesBoundary(*eIt, isIt, localVertCorner) )
                        {
                            std::shared_ptr<BoundaryInteractionVolume> intermediateVolume = std::make_shared<BoundaryInteractionVolume>();
                            BoundaryManager::fillInteractionVolume(intermediateVolume, *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx, /*dummy*/0);
                            intermediateInteractionVolumes_.push_back(intermediateVolume);
                            intermediateVolumesMap_.insert( std::pair<int, int> (globalVertIdx, intermediateVolumesCounter) );
                            intermediateVolumesCounter++;
                        }
                        // otherwise create L-interaction volume
                        else
                        {
                            std::shared_ptr<InteractionVolume> iaVolume1 = std::make_shared<InteractionVolume>(false);
                            std::shared_ptr<InteractionVolume> iaVolume2 = std::make_shared<InteractionVolume>(true);

                            InteractionVolumeManager::fillInteractionVolume(iaVolume1, *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx);
                            InteractionVolumeManager::fillInteractionVolume(iaVolume2, *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx);

                            int faceIndexIn1 = iaVolume1->getFluxFaceIdx();
                            int faceIndexIn2 = iaVolume2->getFluxFaceIdx();

                            typename InteractionVolume::TransmissivityMatrix &T1 = iaVolume1->getTransMatrix();
                            typename InteractionVolume::TransmissivityMatrix &T2 = iaVolume2->getTransMatrix();

                            Scalar t1_0 = T1[faceIndexIn1][0];
                            Scalar t2_0 = T2[faceIndexIn2][0];

                            if ( fabs(t1_0) < fabs(t2_0) )
                                innerInteractionVolumes_.push_back(iaVolume1);
                            else
                                innerInteractionVolumes_.push_back(iaVolume2);

                            std::pair<int, int> iSectionAndNode = std::make_pair(isIt->indexInInside(), localVertCorner);
                            std::pair <int, std::pair<int, int>> elemToPair = std::make_pair(globalElemIdx, iSectionAndNode);
                            innerVolumesMap_.insert(std::pair< std::pair<int, std::pair<int, int> >, int> (elemToPair, innerVolumesCounter));

                            // insert mapping for outer element/facet pair to same interactionvolume
                            // find local vertex index within outside element
                            int localIndexOutside = -1; bool found = false;
                            const ReferenceElement& referenceElementOutside = ReferenceElements::general(isIt->outside().geometry().type());
                            for(int facetCorner = 0; facetCorner < isIt->geometry().corners(); facetCorner++)
                            {
                                // get node index w.r.t outside element
                                int localCornerOutside = referenceElementOutside.subEntity(isIt->indexInOutside(), 1, facetCorner, dim);
                                int globalCornerIdxOutside = problem_().vertexMapper().index( isIt->outside().template subEntity<dim>(localCornerOutside) );
                                if (globalCornerIdxOutside == globalVertIdx)
                                {
                                    localIndexOutside = localCornerOutside;
                                    found = true; break;
                                }
                            }
                            if (found == false)
                                DUNE_THROW(Dune::InvalidStateException, "Couldn't find outside element information!");
                            int outerElementIndex = problem_().elementMapper().index( isIt->outside() );

                            std::pair<int, int> outerIsAndCorner = std::make_pair(isIt->indexInOutside(), localIndexOutside);
                            std::pair <int, std::pair<int, int>> outerElemPair = std::make_pair(outerElementIndex, outerIsAndCorner);
                            innerVolumesMap_.insert(std::pair< std::pair<int, std::pair<int, int> >, int> (outerElemPair, innerVolumesCounter));

                            innerVolumesCounter++;
                        }
                    }
                }
            }
        }

        // Now we handle the remaining elements, i.e. create the remaining MPFA-L volumes
        eIt = problem.gridView().template begin<0>();
        eEndIt = problem.gridView().template end<0>();
        for(; eIt != eEndIt; ++eIt)
        {
            if (eIt->hasBoundaryIntersections())
                continue;
            // global element index
            int globalElemIdx = problem_().elementMapper().index(*eIt);

            // Reference element will be needed later on
            const ReferenceElement& referenceElement = ReferenceElements::general(eIt->geometry().type());

            IntersectionIterator isIt = problem.gridView().ibegin(*eIt);
            IntersectionIterator isEnd = problem.gridView().iend(*eIt);
            for(; isIt != isEnd; ++isIt)
            {
                // loop over the nodes of the facet
                for(int facetNode = 0; facetNode < isIt->geometry().corners(); facetNode++)
                {
                    // get node index w.r.t element
                    int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, facetNode, dim);

                    // get global index
                    int globalVertIdx = problem_().vertexMapper().index( eIt->template subEntity<dim>(localVertCorner) );

                    // only make volume if does not already exist
                    if (!hasVolumeBeenStored(*eIt, isIt->indexInInside(), localVertCorner) )
                    {
                        std::shared_ptr<InteractionVolume> iaVolume1 = std::make_shared<InteractionVolume>(false);
                        std::shared_ptr<InteractionVolume> iaVolume2 = std::make_shared<InteractionVolume>(true);

                        InteractionVolumeManager::fillInteractionVolume(iaVolume1, *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx);
                        InteractionVolumeManager::fillInteractionVolume(iaVolume2, *eIt, problem, isIt->indexInInside(), localVertCorner, globalVertIdx);

                        int FaceIndexin1 = iaVolume1->getFluxFaceIdx();
                        int FaceIndexin2 = iaVolume2->getFluxFaceIdx();

                        typename InteractionVolume::TransmissivityMatrix &T1 = iaVolume1->getTransMatrix();
                        typename InteractionVolume::TransmissivityMatrix &T2 = iaVolume2->getTransMatrix();

                        Scalar t1_0 = T1[FaceIndexin1][0];
                        Scalar t2_0 = T2[FaceIndexin2][0];

                        if ( fabs(t1_0) < fabs(t2_0) )
                            innerInteractionVolumes_.push_back(iaVolume1);
                        else
                            innerInteractionVolumes_.push_back(iaVolume2);

                        std::pair<int, int> iSectionAndNode = std::make_pair(isIt->indexInInside(), localVertCorner);
                        std::pair <int, std::pair<int, int>> elemToPair = std::make_pair(globalElemIdx, iSectionAndNode);
                        innerVolumesMap_.insert(std::pair< std::pair<int, std::pair<int, int> >, int> (elemToPair, innerVolumesCounter));

                        // insert mapping for outer element/facet pair to same interactionvolume
                        // find local vertex index within outside element
                        int localIndexOutside = -1; bool found = false;
                        const ReferenceElement& referenceElementOutside = ReferenceElements::general(isIt->outside().geometry().type());
                        for(int facetCorner = 0; facetCorner < isIt->geometry().corners(); facetCorner++)
                        {
                            // get node index w.r.t element
                            int localCornerOutside = referenceElementOutside.subEntity(isIt->indexInOutside(), 1, facetCorner, dim);
                            int globalCornerIdxOutside = problem_().vertexMapper().index( isIt->outside().template subEntity<dim>(localCornerOutside) );
                            if (globalCornerIdxOutside == globalVertIdx)
                            {
                                localIndexOutside = localCornerOutside;
                                found = true; break;
                            }
                        }
                        if (found == false)
                            DUNE_THROW(Dune::InvalidStateException, "Couldn't find outside element information!");
                        int outerElementIndex = problem_().elementMapper().index( isIt->outside() );

                        std::pair<int, int> outerIsAndCorner = std::make_pair(isIt->indexInOutside(), localIndexOutside);
                        std::pair <int, std::pair<int, int>> outerElemPair = std::make_pair(outerElementIndex, outerIsAndCorner);
                        innerVolumesMap_.insert(std::pair< std::pair<int, std::pair<int, int> >, int> (outerElemPair, innerVolumesCounter));

                        innerVolumesCounter++;
                    }
                }
            }
        }
    }

    // Check if node belongs to a facet which also connects to the boundary
    bool fullFaceTouchesBoundary(const Element &element, IntersectionIterator &isIt, int localNodeIdx)
    {
        // the reference element
        const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

        // the facet itself, thus the node, is not on the boundary. This has been checked in routine "nodeTouchesBoundary"
        // find all the edges that contain this node and
        // check if any of the other nodes of the edges belongs to a boundary facet
        int noEdges = element.subEntities(dim-1);

        for (int edge = 0; edge < noEdges; edge++)
        {
            // an edge ALWAYS contains two nodes
            int otherEdgeNode; bool nodeOnEdge = false;
            for (int edgeCorner = 0; edgeCorner < 2; edgeCorner++)
            {
                int localIndex = referenceElement.subEntity(edge, dim-1, edgeCorner, dim);
                if (localIndex == localNodeIdx)
                {
                    nodeOnEdge = true;
                    if (edgeCorner == 0)
                        otherEdgeNode = referenceElement.subEntity(edge, dim-1, 1, dim);
                    else
                        otherEdgeNode = referenceElement.subEntity(edge, dim-1, 0, dim);
                }
            }

            if (nodeOnEdge == true)
            {
                IntersectionIterator isIt2 = problem_().gridView().ibegin(element);
                IntersectionIterator isIt2End = problem_().gridView().iend(element);
                for(; isIt2 != isIt2End; ++isIt2)
                {
                    if (isIt2 == isIt || !isIt2->boundary())
                        continue;
                    else
                    {
                        for(int corner = 0; corner < isIt2->geometry().corners(); corner++)
                        {
                            // get node index w.r.t element
                            int localVertCorner = referenceElement.subEntity(isIt2->indexInInside(), 1, corner, dim);

                            // if this is true, we found the intersection which also incorporates the node
                            // this is at the same time a boundary facet, as checked by the above if statement
                            if (localVertCorner == otherEdgeNode)
                                return true;
                        }
                    }
                }
            }
        }

        return false;
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
                if (isIt2 == isIt || !isIt2->boundary() )
                    continue;
                else
                    for(int corner = 0; corner < isIt2->geometry().corners(); corner++)
                    {
                        // get node index w.r.t element
                        int localVertCorner = referenceElement.subEntity(isIt2->indexInInside(), 1, corner, dim);

                        if (localVertCorner == localNodeIdx)
                            return true;
                    }
            }
        }

        return false;
    }

    int isBoundaryVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // check the container of the boundary interaction volumes
        // get global vertex index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        std::map<int, int>::const_iterator boundaryIt = boundaryVolumesMap_.find(globalVertIdx);
        if (boundaryIt != boundaryVolumesMap_.end())
            return BoundaryLayers::boundary;
        else
        {
            std::map<int, int>::const_iterator itInter = intermediateVolumesMap_.find(globalVertIdx);
            if( itInter != intermediateVolumesMap_.end())
                return BoundaryLayers::intermediate;
            else
                return BoundaryLayers::interior;
        }
    }

    // Returns a mpfa-L interaction volume from inside the domain
    std::shared_ptr<InteractionVolume> getInnerInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // first check if a corresponding l-Volume exists
        // global Index of the element
        int globalElemIdx = problem_().elementMapper().index( element );

        std::pair<int, int> facetAndVertex = std::make_pair(facetIdxOnElement, localVertIdx);
        std::pair<int, std::pair<int, int> > key = std::make_pair(globalElemIdx, facetAndVertex);

        std::map<std::pair<int, std::pair<int,int> >, int>::const_iterator it = innerVolumesMap_.find(key);
        if (it != innerVolumesMap_.end())
            return innerInteractionVolumes_[it->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Interaction volume not found for " << globalElemIdx << " - " << facetIdxOnElement << " - " << localVertIdx);
    }

    // Returns an mpfa-o interaction volume for a node on the boundary and an equation index
    std::shared_ptr<BoundaryInteractionVolume> getBoundaryInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx, int phaseIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        // first check the boundary interaction volumes, then the intermediate volumes
        std::map<int, int>::const_iterator itBoundary = boundaryVolumesMap_.find(globalVertIdx);
        if (itBoundary != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[itBoundary->second][phaseIdx];
        else
           DUNE_THROW(Dune::InvalidStateException, "Boundary interaction volume not found!!!");
    }
    std::shared_ptr<BoundaryInteractionVolume> getBoundaryInteractionVolumeFromGlobalVertexIndex(int globalVertIdx, int phaseIdx) const
    {
        std::map<int, int>::const_iterator itBoundary = boundaryVolumesMap_.find(globalVertIdx);
        if (itBoundary != boundaryVolumesMap_.end())
            return boundaryInteractionVolumes_[itBoundary->second][phaseIdx];
        else
            DUNE_THROW(Dune::InvalidStateException, "Boundary interaction volume not found!!!");
    }

    // Function for returning the intermediate (o-type) interaction volumes
    std::shared_ptr<BoundaryInteractionVolume> getInsideInteractionVolume(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // get global node index
        int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

        std::map<int, int>::const_iterator itInter = intermediateVolumesMap_.find(globalVertIdx);
        if( itInter != intermediateVolumesMap_.end())
            return intermediateInteractionVolumes_[itInter->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Intermediate interaction volume not found!!!");
    }
    std::shared_ptr<BoundaryInteractionVolume> getInsideInteractionVolumeFromGlobalVertexIndex(int globalVertIdx, int phaseIdx) const
    {
        std::map<int, int>::const_iterator itInter = intermediateVolumesMap_.find(globalVertIdx);
        if( itInter != intermediateVolumesMap_.end())
            return intermediateInteractionVolumes_[itInter->second];
        else
            DUNE_THROW(Dune::InvalidStateException, "Intermediate interaction volume not found!!!");
    }


    bool hasVolumeBeenStored(const Element &element, int facetIdxOnElement, int localVertIdx) const
    {
        // first check the container of the inner interaction volumes
        int globalElemIdx = problem_().elementMapper().index(element);

        std::pair<int, int> facetAndVertex = std::make_pair(facetIdxOnElement, localVertIdx);
        std::pair<int, std::pair<int,int> > key = std::make_pair(globalElemIdx, facetAndVertex);

        std::map< std::pair<int, std::pair<int, int> >, int >::const_iterator it = innerVolumesMap_.find(key);

        if (it != innerVolumesMap_.end())
            return true;
        else
        {
            // check the container of the boundary volumes
            // get global vertex index
            int globalVertIdx = problem_().vertexMapper().index( element.template subEntity<dim>(localVertIdx) );

            std::map<int, int>::const_iterator itBoundary = boundaryVolumesMap_.find(globalVertIdx);
            if (itBoundary != boundaryVolumesMap_.end())
                return true;
            else
            {
                std::map<int, int>::const_iterator itInter = intermediateVolumesMap_.find(globalVertIdx);
                if( itInter != intermediateVolumesMap_.end())
                    return true;
                else
                    return false;
            }
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
    std::vector< std::shared_ptr<BoundaryInteractionVolume> > intermediateInteractionVolumes_;

    std::map< int, int > boundaryVolumesMap_;
    std::map< std::pair< int, std::pair<int, int> >, int> innerVolumesMap_;
    std::map< int, int > intermediateVolumesMap_;
};
}
#endif /* DUMUX_MPFA_O_INTERACTIONVOLUMECONTAINER_HH */
