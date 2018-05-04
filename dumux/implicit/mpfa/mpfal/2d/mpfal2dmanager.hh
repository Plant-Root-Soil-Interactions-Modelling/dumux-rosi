#ifndef DUMUX_MPFALINTERACTIONVOLUMEMANAGER2D_HH
#define DUMUX_MPFALINTERACTIONVOLUMEMANAGER2D_HH

#include <dumux/common/math.hh>

namespace Dumux{

template<class TypeTag>
class MpfaL2DManager
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryLayers) BoundaryLayers;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;


public:

    static void fillInteractionVolume(std::shared_ptr<InteractionVolume> interactionVolume,
                                            const Element &element, const Problem &problem,
                                                int facetIdxOnElement, int nodeIndex, int centralVertGlobalIdx)
    {
        if (dim != 2)
            DUNE_THROW(Dune::NotImplemented, "MPFA-O 2D interactionvolumefiller has been called for dim = " << dim);

        // pass the vertex & element mapper
        const VertexMapper &vertexMapper = problem.vertexMapper();
        const ElementMapper &elementMapper = problem.elementMapper();

        // introduce matrix R for vector rotation and R is initialized as zero matrix
        Dune::FieldMatrix<Scalar, dim, dim> R(0);

        // a vector initialized with zeros, used to switch sign of vector: v -> -v
        const GlobalPosition helpVector(0);

        // Fill rotation matrix R
        R[0][1] = 1;
        R[1][0] = -1;

        // get information on the actual element
        const GlobalPosition& Element1Position = element.geometry().center();
        const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

        // interaction volume is created around a vertex of the mesh, get coordinates
        const GlobalPosition& centralVertPosition = element.geometry().corner(nodeIndex);

        // find the two intersections of the actual element that contain this node
        std::vector<int> IsIndicesInside(2, -1);
        int counter = 0;
        IntersectionIterator isIt = problem.gridView().ibegin(element);
        IntersectionIterator isItEnd = problem.gridView().iend(element);
        // declaration of the two intersectioniterators
        // later the respective intersections will be assigned to them
        // errors will be caught below in case something goes wrong
        IntersectionIterator isIt12temp = isIt;
        IntersectionIterator isIt14temp = isIt;
        for (; isIt != isItEnd; ++isIt)
        {
            if(counter == 2)
                break;

            // get the global coordinate and global vertex index of corners
            for (int i = 0; i < isIt->geometry().corners(); ++i)
            {
                int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, i, dim);
                int globalVertIdxcorner = vertexMapper.index( element.template subEntity<dim>(localVertCorner) );

                if (globalVertIdxcorner == centralVertGlobalIdx && counter == 0)
                {
                    IsIndicesInside[counter] = isIt->indexInInside();
                    isIt12temp = isIt;
                    counter++;
                    break;
                }

                if (globalVertIdxcorner == centralVertGlobalIdx && counter == 1)
                {
                    IsIndicesInside[counter] = isIt->indexInInside();
                    isIt14temp = isIt;
                    counter++;
                    break;
                }
            }
        }

        if (counter != 2)
            DUNE_THROW(Dune::NotImplemented, "Error! counter != 2");
        if(IsIndicesInside[0] == -1 || IsIndicesInside[1] == -1)
            DUNE_THROW(Dune::NotImplemented, "Error! Facets not found!");
        if(IsIndicesInside[0] != facetIdxOnElement && IsIndicesInside[1] != facetIdxOnElement)
            DUNE_THROW(Dune::NotImplemented, "Facet doesn't correspond to passed facet index");

        // now we check in which order the center and the two continuity points
        // on the intersections form a right handed system
        // the connecting vectors between element center and continuity points are stored
        std::vector<GlobalPosition> connectVectors(2);
        std::vector<GlobalPosition> contiPoints(2);
        bool swap = calcConnectionVectors(element, IsIndicesInside, centralVertPosition, connectVectors, contiPoints);

        // Now we have the two intersections of the first element
        // in case the two intersections were swapped in order to
        // get a right handed system we have to swap the two intersections here, too
        IntersectionIterator isIt12 = isIt12temp;
        IntersectionIterator isIt14 = isIt14temp;
        if (swap)
        {
            isIt12 = isIt14temp;
            isIt14 = isIt12temp;
        }

        if (!interactionVolume->secondTriangle())
        {
            // set the flag indicating an interactionvolume has been (or rather is about to be) stored
            interactionVolume->setStored();

            if (isIt12->indexInInside() == facetIdxOnElement)
                interactionVolume->setFluxFaceIndex(0);
            else if (isIt14->indexInInside() == facetIdxOnElement)
                interactionVolume->setFluxFaceIndex(1);
            else
                DUNE_THROW(Dune::NotImplemented, "None of the Facets correspond to passed facet index");

            // compute normal vectors nu1 & nu2
            GlobalPosition nu2(0);
            R.mv(helpVector - connectVectors[0], nu2);
            interactionVolume->setNu(nu2, 1);

            GlobalPosition nu1(0);
            R.mv(connectVectors[1], nu1);
            interactionVolume->setNu(nu1, 0);

            GlobalPosition nu7(0);
            R.mv(centralVertPosition - Element1Position, nu7);
            interactionVolume->setNu(nu7, 6);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition Rnu1(0);
            R.umv(nu1, Rnu1);
            Scalar T1 = fabs(nu2 * Rnu1);
            interactionVolume->setT(T1, 0);
            interactionVolume->setSubVolume(element, 0);

            GlobalPosition normal1 = isIt12->centerUnitOuterNormal();
            normal1 *= isIt12->geometry().volume() / 2;
            interactionVolume->setNormal(normal1, 0);
            interactionVolume->setFaceToFacetMaps(0, isIt12->indexInInside(), 0);
            interactionVolume->setFaceToFacetMaps(1, isIt12->indexInOutside(), 0);

            GlobalPosition normal2 = isIt14->centerUnitOuterNormal();
            normal2 *= isIt14->geometry().volume() / 2;
            interactionVolume->setNormal(normal2, 1);
            interactionVolume->setFaceToFacetMaps(0, isIt14->indexInInside(), 1);
            interactionVolume->setFaceToFacetMaps(2, isIt14->indexInOutside(), 1);

            // get second element and its center
            interactionVolume->setSubVolume(isIt12->outside(), 1);
            const GlobalPosition element2Center = isIt12->outside().geometry().center();

            GlobalPosition X3 = contiPoints[0] - element2Center;
            GlobalPosition nu3(0);
            R.mv(X3, nu3);
            interactionVolume->setNu(nu3, 2);

            GlobalPosition X4 = centralVertPosition - element2Center;
            GlobalPosition nu4(0);
            R.mv(helpVector - X4, nu4);
            interactionVolume->setNu(nu4, 3);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition RX4(0);
            R.umv(X4, RX4);
            Scalar T2 = fabs(X3 * RX4);
            interactionVolume->setT(T2, 1);

            // get third element and its center
            interactionVolume->setSubVolume(isIt14->outside(), 2);
            const GlobalPosition element3Center = isIt14->outside().geometry().center();

            GlobalPosition X5 = centralVertPosition - element3Center;
            GlobalPosition nu5(0);
            R.mv(X5, nu5);
            interactionVolume->setNu(nu5, 4);

            GlobalPosition X6 = contiPoints[1] - element3Center;
            GlobalPosition nu6(0);
            R.mv(helpVector - X6, nu6);
            interactionVolume->setNu(nu6, 5);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition RX5(0);
            R.umv(X5, RX5);
            Scalar T3 = fabs(X6 * RX5);
            interactionVolume->setT(T3, 2);
        }
        else
        {
            IntersectionIterator facetIS = isIt12;
            interactionVolume->setFluxFaceIndex(1);
            if(isIt14->indexInInside() == facetIdxOnElement)
            {
                facetIS = isIt14;
                interactionVolume->setFluxFaceIndex(0);
            }

            // find the two intersections of the element that contain this node
            std::vector<int> IsIndicesInside2(2, -1);

            const GlobalPosition centralElemCenter = facetIS->outside().geometry().center();

            int counter2 = 0;
            IntersectionIterator isItElem2 = problem.gridView().ibegin( facetIS->outside() );
            IntersectionIterator isItElem2End = problem.gridView().iend( facetIS->outside() );

            // declaration of the two intersectioniterators
            // later the respective intersections will be assigned to them
            // errors will be caught below in case something goes wrong
            IntersectionIterator isItElem2temp1 = isItElem2;
            IntersectionIterator isItElem2temp2 = isItElem2;
            for (; isItElem2 != isItElem2End; ++isItElem2)
            {
                if(counter2 == 2)
                    break;

                // get the global coordinate and global vertex index of corner1234
                for (int i = 0; i < isItElem2->geometry().corners(); ++i)
                {

                    int localVertCorner = referenceElement.subEntity(isItElem2->indexInInside(), 1, i, dim);
                    int globalVertIdxcorner = vertexMapper.index( facetIS->outside()->template subEntity<dim>(localVertCorner) );

                    if (globalVertIdxcorner == centralVertGlobalIdx && counter2 == 0)
                    {
                        IsIndicesInside2[counter2] = isItElem2->indexInInside();
                        isItElem2temp1 = isItElem2;
                        counter2++;
                        break;
                    }

                    if (globalVertIdxcorner == centralVertGlobalIdx && counter2 == 1)
                    {
                        IsIndicesInside2[counter2] = isItElem2->indexInInside();
                        isItElem2temp2 = isItElem2;
                        counter2++;
                        break;
                    }
                }
            }

            if (counter2 != 2)
                DUNE_THROW(Dune::NotImplemented, "Error! counter2 != 2");
            if(IsIndicesInside2[0] == -1 || IsIndicesInside2[1] == -1)
                DUNE_THROW(Dune::NotImplemented, "Error! Facets not found!");
            if(IsIndicesInside2[0] != facetIS->indexInOutside() && IsIndicesInside2[1] != facetIS->indexInOutside())
                DUNE_THROW(Dune::NotImplemented, "Facet doesn't correspond to passed facet index");

            // now we check in which order the center of the element and the two continuity points
            // on the intersections form a right handed system
            // the connecting vectors between element center and continuity points are stored
            std::vector<GlobalPosition> connectVectors2(2);
            std::vector<GlobalPosition> contiPoints2(2);
            bool swap = calcConnectionVectors(facetIS->outside(), IsIndicesInside2, centralVertPosition, connectVectors2, contiPoints2);

            // Now we have the two intersections of the first element
            // in case the two intersections were swapped in order to
            // get a right handed system we have to swap the two intersections here, too
            IntersectionIterator isItElem2_1 = isItElem2temp1;
            IntersectionIterator isItElem2_2 = isItElem2temp2;
            if (swap)
            {
                isItElem2_1 = isItElem2temp2;
                isItElem2_2 = isItElem2temp1;
            }

            // set the flag indicating an interactionvolume has been (or rather is about to be) stored
            interactionVolume->setStored();

            // compute normal vectors nu1 & nu2
            GlobalPosition nu2(0);
            R.mv(helpVector - connectVectors2[0], nu2);
            interactionVolume->setNu(nu2, 1);

            GlobalPosition nu1(0);
            R.mv(connectVectors2[1], nu1);
            interactionVolume->setNu(nu1, 0);

            GlobalPosition nu7(0);
            R.mv(centralVertPosition - centralElemCenter, nu7);
            interactionVolume->setNu(nu7, 6);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition Rnu1(0);
            R.umv(nu1, Rnu1);
            Scalar T1 = fabs(nu2 * Rnu1);
            interactionVolume->setT(T1, 0);
            interactionVolume->setSubVolume(facetIS->outside(), 0);

            GlobalPosition normal1 = isItElem2_1->centerUnitOuterNormal();
            normal1 *= isItElem2_1->geometry().volume() / 2;
            interactionVolume->setNormal(normal1, 0);
            interactionVolume->setFaceToFacetMaps(0, isItElem2_1->indexInInside(), 0);
            interactionVolume->setFaceToFacetMaps(1, isItElem2_1->indexInOutside(), 0);

            GlobalPosition normal2 = isItElem2_2->centerUnitOuterNormal();
            normal2 *= isItElem2_2->geometry().volume() / 2;
            interactionVolume->setNormal(normal2, 1);
            interactionVolume->setFaceToFacetMaps(0, isItElem2_2->indexInInside(), 1);
            interactionVolume->setFaceToFacetMaps(2, isItElem2_2->indexInOutside(), 1);

            // get second element and its center
            interactionVolume->setSubVolume(isItElem2_1->outside(), 1);
            const GlobalPosition element2Center = isItElem2_1->outside().geometry().center();

            GlobalPosition X3 = contiPoints2[0] - element2Center;
            GlobalPosition nu3(0);
            R.mv(X3, nu3);
            interactionVolume->setNu(nu3, 2);

            GlobalPosition X4 = centralVertPosition - element2Center;
            GlobalPosition nu4(0);
            R.mv(helpVector - X4, nu4);
            interactionVolume->setNu(nu4, 3);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition RX4(0);
            R.umv(X4, RX4);
            Scalar T2 = fabs(X3 * RX4);
            interactionVolume->setT(T2, 1);

            // get third element and its center
            interactionVolume->setSubVolume(isItElem2_2->outside(), 2);
            const GlobalPosition element3Center = isItElem2_2->outside().geometry().center();

            GlobalPosition X5 = centralVertPosition - element3Center;
            GlobalPosition nu5(0);
            R.mv(X5, nu5);
            interactionVolume->setNu(nu5, 4);

            GlobalPosition X6 = contiPoints2[1] - element3Center;
            GlobalPosition nu6(0);
            R.mv(helpVector - X6, nu6);
            interactionVolume->setNu(nu6, 5);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition RX6(0);
            R.umv(X6, RX6);
            Scalar T3 = fabs(X5 * RX6);
            interactionVolume->setT(T3, 2);
        }
        // uncomment line below to have the interactionvolume printed out
        //interactionVolume->printInteractionVolume();
        // set up transmissivity matrix and store in interaction volume
        interactionVolume->calculateTransmissivityMatrix(problem);
    } // end of function


    // fills a storage vector with the elements of the stencil
    // assumes that the first element of the input container
    // is the element of which the stencil is desired
    static int getElementsOfStencil(std::vector<Element> &container, const Problem &problem)
    {
        int numNeighbors = 1;

        // Reference element will be needed later on
        const ReferenceElement& referenceElement = ReferenceElements::general(container[0].geometry().type());

        // Loop over the facets and nodes of the element and the corresponding interaction volumes
        // to fill the the container with info about the elements of the stencil
        IntersectionIterator isIt = problem.gridView().ibegin(container[0]);
        IntersectionIterator isEnd = problem.gridView().iend(container[0]);
        for(; isIt != isEnd; ++isIt)
        {
            for (int nodeIdx = 0; nodeIdx < isIt->geometry().corners(); nodeIdx++)
            {
                int localVertIdx = referenceElement.subEntity(isIt->indexInInside(), 1, nodeIdx, dim);

                std::vector<Element> tmpElements;
                int boundaryLayer = problem.model().interactionVolumeContainer().isBoundaryVolume(*(container[0]), isIt->indexInInside(), localVertIdx);
                if(boundaryLayer == BoundaryLayers::interior)
                    problem.model().getInnerInteractionVolume(container[0], isIt->indexInInside(), localVertIdx)->passElementsInRegion(tmpElements);
                else if (boundaryLayer == BoundaryLayers::intermediate)
                    problem.model().getInsideInteractionVolume(container[0], isIt->indexInInside(), localVertIdx)->passElementsInRegion(tmpElements);
                else
                    problem.model().getBoundaryInteractionVolume(container[0], isIt->indexInInside(), localVertIdx, /*dummy*/0)->passElementsInRegion(tmpElements);

                // insert new elements in provided container
                for (int newElIdx = 0; newElIdx < tmpElements.size(); newElIdx++)
                {
                    // first check if element has already been stored
                    bool found = false;
                    for (int neighIdx = 0; neighIdx < numNeighbors; neighIdx++)
                    {
                        if (tmpElements[newElIdx] == container[neighIdx])
                        {
                            found = true; break;
                        }
                    }
                    // insert new element only if found = false
                    if (found == false)
                    {
                        container.push_back(tmpElements[newElIdx]);
                        numNeighbors++;
                    }
                }
            }
        }

        return numNeighbors;
    }


    static bool calcConnectionVectors(const Element& element,
                                      std::vector<int>& IsIndicesInsideElement,
                                      const GlobalPosition corner,
                                      std::vector<GlobalPosition>& connectVectorStorage,
                                      std::vector<GlobalPosition>& contiPointStorage)
    {
        bool swap = false;

        Scalar continuityPoint = GET_PROP_VALUE(TypeTag, MpfaContinuityPoint);
        const GlobalPosition isCenter1 = element.template subEntity<1>(IsIndicesInsideElement[0])->geometry().center();
        const GlobalPosition isCenter2 = element.template subEntity<1>(IsIndicesInsideElement[1])->geometry().center();

        // calculate the positions on which pressure continuity is demanded
        GlobalPosition contiPoint1 = isCenter1;
        GlobalPosition contiPoint2 = isCenter2;

        // calculate connection vectors between element center and continuity points
        const GlobalPosition elCenter = element.geometry().center();
        GlobalPosition X1 = contiPoint1;
        X1 -= elCenter;
        GlobalPosition X2 = contiPoint2;
        X2 -= elCenter;

        Scalar product = X1[0]*X2[1] - X1[1]*X2[0];

        // if vector product < 0, we don't have a right handed system and we have to swap the intersection indices
        if (product < 0)
        {
            swap = true;
            int swapIndex = IsIndicesInsideElement[0];
            IsIndicesInsideElement[0] = IsIndicesInsideElement[1];
            IsIndicesInsideElement[1] = swapIndex;

            connectVectorStorage[0] = X2;
            connectVectorStorage[1] = X1;
            contiPointStorage[0] = contiPoint2;
            contiPointStorage[1] = contiPoint1;
        }
        else if (product > 0)
        {
            connectVectorStorage[0] = X1;
            connectVectorStorage[1] = X2;
            contiPointStorage[0] = contiPoint1;
            contiPointStorage[1] = contiPoint2;
        }
        else
            DUNE_THROW(Dune::MathError, "product = 0!");

        return swap;
    }
};
}
#endif
