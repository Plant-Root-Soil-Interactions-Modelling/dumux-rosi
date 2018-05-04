#ifndef DUMUX_MPFAOINTERACTIONVOLUMEMANAGER2D_HH
#define DUMUX_MPFAOINTERACTIONVOLUMEMANAGER2D_HH

#include <dumux/common/math.hh>

namespace Dumux{

template<class TypeTag>
class MpfaO2DManager
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP(TypeTag, MpfaMethods) MpfaMethods;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryLayers) BoundaryLayers;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaFaceTypes) FaceTypes;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        pressureIdx = Indices::pressureIdx
    };
    typedef Dumux::MpfaO2DInteractionVolume<TypeTag> InteractionVolume;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename InteractionVolume::SubVolume SubVolume;
    typedef typename InteractionVolume::SubVolumeFace SubVolumeFace;

public:

    static void fillInteractionVolume(std::shared_ptr<InteractionVolume> interactionVolume,
                                            const Element &element, const Problem &problem,
                                                int facetIdxOnElement, int nodeIndex, int centralVertGlobalIdx, int phaseIdx)
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

        // evaluate matrix R
        R[0][1] = 1;
        R[1][0] = -1;

        // get global Index of cell 1
        int element1GlobalIdx = elementMapper.index(element);

        // Reference element will be needed later on
        const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

        // Build interaction volume for the node of the actual element
        // boolean to indicate whether interaction volume creation has been finished
        bool volumeCreated = false;

        // counters for the elements/faces that have been added to the interaction region.
        int elementCounter = 0; int faceCounter = 0;

        // get coordinate of the actual node
        const GlobalPosition& centralVertPosition = element.geometry().corner(nodeIndex);

        // find the two intersections of the element that contain this node
        std::vector<int> IsIndicesInside1(2, -1);

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

            // get the global coordinate and global vertex index of corner1234
            for (int i = 0; i < isIt->geometry().corners(); ++i)
            {

                int localVertCorner = referenceElement.subEntity(isIt->indexInInside(), 1, i, dim);
                int globalVertIdxcorner = vertexMapper.index( element.template subEntity<dim>(localVertCorner) );

                if (globalVertIdxcorner == centralVertGlobalIdx && counter == 0)
                {
                    IsIndicesInside1[counter] = isIt->indexInInside();
                    isIt12temp = isIt;
                    counter++;
                    break;
                }

                if (globalVertIdxcorner == centralVertGlobalIdx && counter == 1)
                {
                    IsIndicesInside1[counter] = isIt->indexInInside();
                    isIt14temp = isIt;
                    counter++;
                    break;
                }
            }
        }

        if (counter != 2)
            DUNE_THROW(Dune::InvalidStateException, "Error in mpfaointeractionvolume2dfiller.hh, l.132");
        if(IsIndicesInside1[0] == -1 || IsIndicesInside1[1] == -1)
            DUNE_THROW(Dune::InvalidStateException, "Error in mpfaointeractionvolume2dfiller.hh, l.134");

        // now we check in which order the center and the two continuity points
        // on the intersections form a right handed system
        // the connecting vectors between element center and continuity points are stored
        std::vector<GlobalPosition> connectVectors1(2);
        std::vector<GlobalPosition> contiPoints1(2);
        bool swap = calcConnectionVectors(element, IsIndicesInside1, centralVertPosition, connectVectors1, contiPoints1);

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

        // set the flag indicating an interactionvolume has been (or rather is about to be) stored
        interactionVolume->setStored();
        interactionVolume->setCentralVertex(centralVertPosition);
        interactionVolume->setCentralVertexIndex(centralVertGlobalIdx);

        // compute normal vectors nu11,nu12
        GlobalPosition nu12(0);
        R.mv(helpVector - connectVectors1[0], nu12);

        GlobalPosition nu11(0);
        R.mv(connectVectors1[1], nu11);

        // compute T, the area of quadrilateral made by normal vectors 'nu'
        GlobalPosition Rnu11(0);
        R.umv(nu11, Rnu11);
        Scalar T1 = fabs(nu12 * Rnu11);

        // Make first subVolume object
        SubVolume subVol1;
        subVol1.element = element;
        subVol1.globalIdx = element1GlobalIdx;
        subVol1.localIsIndex[0] = IsIndicesInside1[0];
        subVol1.localIsIndex[1] = IsIndicesInside1[1];
        subVol1.contiPoints[0] = contiPoints1[0];
        subVol1.contiPoints[1] = contiPoints1[1];
        subVol1.nu[0] = nu11;
        subVol1.nu[1] = nu12;
        subVol1.T = T1;

        // pass sub volume to interactionVolume container and increment counter
        interactionVolume->addSubVolume(subVol1);
        interactionVolume->setElementToSubVolMap(element1GlobalIdx,elementCounter);
        elementCounter++;

        IntersectionIterator nextIS = isIt12;

        // we rotate counter clockwise and create new sub volumes/faces
        while (!nextIS->boundary())
        {
            // get volume of face
            IntersectionGeometry isGeometry = nextIS->geometry();
            Scalar faceVol = isGeometry.volume() / 2.0;
            // get outer normal vector of face
            GlobalPosition unitOuterNormal = nextIS->centerUnitOuterNormal();
            // get center of the intersection
            GlobalPosition center = isGeometry.center();

            // check if we are on an interior Boundary
            const int faceType = problem.isInternalBoundary(*nextIS);

            // Check if we are on the interface to the first element
            if (nextIS->outside() == element)
            {
                // Make the last sub volume face object
                SubVolumeFace subFace;
                subFace.normal = unitOuterNormal;
                subFace.center = center;
                subFace.area = faceVol;
                subFace.positiveSubVolIndex = elementCounter - 1;
                subFace.negativeSubVolIndex = 0;
                subFace.localIdxOnPositive = 0;
                subFace.localIdxOnNegative = 1;
                subFace.onBoundary = false;
                subFace.faceType = faceType;
                // if internal dirichlet face, add face without face unknown
                // otherwise, add face with face unknown
                if (faceType == FaceTypes::InternalDirichletFace)
                    interactionVolume->addDirichletFace(faceCounter);
                else
                    interactionVolume->addInteriorFace(faceCounter);

                // pass sub volume face to interactionVolume container
                interactionVolume->addSubVolumeFace(subFace);

                // also set the mapping from the local face indices of the sub volumes to the face
                interactionVolume->setSubVolToFaceMap(subFace.positiveSubVolIndex, subFace.localIdxOnPositive, faceCounter);

                if (faceType != FaceTypes::InternalDirichletFace && faceType != FaceTypes::InternalFluxFace)
                    interactionVolume->setSubVolToFaceMap(subFace.negativeSubVolIndex, subFace.localIdxOnNegative, faceCounter);

                // if facet coupling is on, move the conti point back to the center (in case it was moved)
                if ( (faceType == FaceTypes::InternalDirichletFace || faceType == FaceTypes::InternalFluxFace) && GET_PROP_VALUE(TypeTag, FacetCoupling) )
                {
                    GlobalPosition newContiPointPositive = center;

                    // get connection vector from element to conti point
                    GlobalPosition newConnectVec = newContiPointPositive;
                    newConnectVec -= interactionVolume->getSubVolume(subFace.positiveSubVolIndex).element.geometry().center();

                    // calculate the new nu2
                    GlobalPosition newNu2(0);
                    R.mv(helpVector - newConnectVec, newNu2);

                    // compute new T, the area of quadrilateral made by normal vectors 'nu'
                    GlobalPosition Rnu1(0);
                    GlobalPosition oldNu1 = interactionVolume->getSubVolume(subFace.positiveSubVolIndex).nu[0];
                    R.umv(oldNu1, Rnu1);
                    Scalar newT = fabs(newNu2 * Rnu1);

                    // update contiPoint for local index 0, nu for local index 1 and T
                    interactionVolume->updateSubVolumeData(subFace.positiveSubVolIndex, newContiPointPositive, 0, newNu2, 1, newT);
                }
                faceCounter++;

                // if we are on an internal Flux or Dirichlet boundary we have to duplicate the face swapping the sub volume indices
                if (faceType == FaceTypes::InternalDirichletFace || faceType == FaceTypes::InternalFluxFace)
                {
                    SubVolumeFace subFaceDouble;
                    subFaceDouble.normal = helpVector - unitOuterNormal;
                    subFaceDouble.center = center;
                    subFaceDouble.area = faceVol;
                    subFaceDouble.positiveSubVolIndex = 0;
                    subFaceDouble.negativeSubVolIndex = elementCounter - 1;
                    subFaceDouble.localIdxOnPositive = 1;
                    subFaceDouble.localIdxOnNegative = 0;
                    subFaceDouble.onBoundary = false;
                    subFaceDouble.faceType = faceType;
                    if (faceType == FaceTypes::InternalFluxFace)
                        interactionVolume->addInteriorFace(faceCounter);
                    else
                        interactionVolume->addDirichletFace(faceCounter);
                    // pass sub volume face to interactionVolume container and increment counter
                    interactionVolume->addSubVolumeFace(subFaceDouble);
                    // also set the mapping from the local face indices of the sub volumes to the face
                    interactionVolume->setSubVolToFaceMap(subFaceDouble.positiveSubVolIndex, subFaceDouble.localIdxOnPositive, faceCounter);

                    // if facet coupling is on, move the conti point of the double face back to the center (in case it was moved)
                    if ( GET_PROP_VALUE(TypeTag, FacetCoupling) )
                    {
                        GlobalPosition newContiPointPositive = center;
                        // get connection vector from element to conti point
                        GlobalPosition newConnectVec = newContiPointPositive;
                        newConnectVec -= interactionVolume->getSubVolume(subFaceDouble.positiveSubVolIndex).element.geometry().center();

                        // calculate the new nu2
                        GlobalPosition newNu1(0);
                        R.mv(newConnectVec, newNu1);

                        // compute new T, the area of quadrilateral made by normal vectors 'nu'
                        GlobalPosition Rnu1(0);
                        GlobalPosition oldNu2 = interactionVolume->getSubVolume(subFaceDouble.positiveSubVolIndex).nu[1];
                        R.umv(newNu1, Rnu1);
                        Scalar newT = fabs(oldNu2 * Rnu1);

                        // update contiPoint for local index 1, nu for local index 0 and T
                        interactionVolume->updateSubVolumeData(subFaceDouble.positiveSubVolIndex, newContiPointPositive, 1, newNu1, 0, newT);
                    }
                    faceCounter++;
                }

                // If we are here, that means the outer element is the first one again, i.e. ia-volume creation has finished
                volumeCreated = true;
                interactionVolume->setInteriorVolume(true);
                break;
            }
            else
            {
                // Make sub volume face object
                SubVolumeFace subFace;
                subFace.normal = unitOuterNormal;
                subFace.center = center;
                subFace.area = faceVol;
                subFace.positiveSubVolIndex = elementCounter - 1;
                subFace.negativeSubVolIndex = elementCounter;
                subFace.localIdxOnPositive = 0;
                subFace.localIdxOnNegative = 1;
                subFace.onBoundary = false;
                subFace.faceType = faceType;
                // if internal dirichlet face, add face without face unknown
                // otherwise, add face with face unknown
                if (faceType == FaceTypes::InternalDirichletFace)
                    interactionVolume->addDirichletFace(faceCounter);
                else
                    interactionVolume->addInteriorFace(faceCounter);

                // pass sub volume face to interactionVolume container and increment counter
                interactionVolume->addSubVolumeFace(subFace);

                interactionVolume->setSubVolToFaceMap(subFace.positiveSubVolIndex, subFace.localIdxOnPositive, faceCounter);
                if (faceType != FaceTypes::InternalDirichletFace && faceType != FaceTypes::InternalFluxFace)
                    interactionVolume->setSubVolToFaceMap(subFace.negativeSubVolIndex, subFace.localIdxOnNegative, faceCounter);

                // if facet coupling is on, move the conti point of the double face back to the center (in case it was moved)
                if ( (faceType == FaceTypes::InternalFluxFace || faceType == FaceTypes::InternalDirichletFace) && GET_PROP_VALUE(TypeTag, FacetCoupling) )
                {
                    GlobalPosition newContiPointPositive = center;

                    // get connection vector from element to conti point
                    GlobalPosition newConnectVec = newContiPointPositive;
                    newConnectVec -= interactionVolume->getSubVolume(subFace.positiveSubVolIndex).element.geometry().center();

                    // calculate the new nu2
                    GlobalPosition newNu2(0);
                    R.mv(helpVector - newConnectVec, newNu2);

                    // compute new T, the area of quadrilateral made by normal vectors 'nu'
                    GlobalPosition Rnu1(0);
                    GlobalPosition oldNu1 = interactionVolume->getSubVolume(subFace.positiveSubVolIndex).nu[0];
                    R.umv(oldNu1, Rnu1);
                    Scalar newT = fabs(newNu2 * Rnu1);

                    // update contiPoint for local index 0, nu for local index 1 and T
                    interactionVolume->updateSubVolumeData(subFace.positiveSubVolIndex, newContiPointPositive, 0, newNu2, 1, newT);
                }
                faceCounter++;

                // if we are on an internal Flux boundary we have to duplicate the face swapping the sub volume indices
                if (faceType == FaceTypes::InternalDirichletFace || faceType == FaceTypes::InternalFluxFace)
                {
                    SubVolumeFace subFaceDouble;
                    subFaceDouble.normal = helpVector - unitOuterNormal;
                    subFaceDouble.center = center;
                    subFaceDouble.area = faceVol;
                    subFaceDouble.positiveSubVolIndex = elementCounter;
                    subFaceDouble.negativeSubVolIndex = elementCounter - 1;
                    subFaceDouble.localIdxOnPositive = 1;
                    subFaceDouble.localIdxOnNegative = 0;
                    subFaceDouble.onBoundary = false;
                    subFaceDouble.faceType = faceType;
                    if (faceType == FaceTypes::InternalFluxFace)
                        interactionVolume->addInteriorFace(faceCounter);
                    else
                        interactionVolume->addDirichletFace(faceCounter);
                    // pass sub volume face to interactionVolume container and increment counter
                    interactionVolume->addSubVolumeFace(subFaceDouble);
                    // also set the mapping from the local face indices of the sub volumes to the face
                    interactionVolume->setSubVolToFaceMap(subFaceDouble.positiveSubVolIndex, subFaceDouble.localIdxOnPositive, faceCounter);
                    faceCounter++;
                }
            }

            // Eventual bug tracking... TODO: to be erased after testing
            if (volumeCreated)
            {
                std::cout << "I'm still here even though i created the volume!?" << std::endl;
                break;
            }

            // obtain the information of the outside element
            Element newElement = nextIS->outside();
            int newElementGlobalIdx = elementMapper.index(newElement);

            // find the intersection of the new element that contains the central node
            IntersectionIterator newIS;
            bool found = false;
            IntersectionIterator isNewEl = problem.gridView().ibegin(newElement);
            for(; isNewEl != problem.gridView().iend(newElement); ++isNewEl)
            {
                if (found == true)
                    break;

                for (int i = 0; i < isNewEl->geometry().corners(); ++i)
                {
                    int localVertCorner = referenceElement.subEntity(isNewEl->indexInInside(), 1, i, dim);

                    int globalVertIdxcorner = vertexMapper.index( newElement.template subEntity<dim>(localVertCorner) );

                    if (globalVertIdxcorner == centralVertGlobalIdx && isNewEl->indexInInside() != nextIS->indexInOutside())
                    {
                        newIS = isNewEl; found = true;
                        break;
                    }
                }
            }
            if (found == false)
                DUNE_THROW(Dune::InvalidStateException, "Couldn't find facet of the neighbouring element that contains the same node!");

            // copy the two intersection indices in a container and calculate
            // continuity points. Boolean swap used as an error indicator
            std::vector<int> IsIndicesInside(2);
            IsIndicesInside[0] = newIS->indexInInside();
            IsIndicesInside[1] = nextIS->indexInOutside();

            std::vector<GlobalPosition> connectVectors(2);
            std::vector<GlobalPosition> contiPoints(2);
            bool swap = calcConnectionVectors(newElement, IsIndicesInside, centralVertPosition, connectVectors, contiPoints);
            if (swap)
                DUNE_THROW(Dune::InvalidStateException, "This swapping shouldn't happen!");

            // compute normal vectors nu1,nu2
            GlobalPosition nu2(0);
            R.mv(helpVector - connectVectors[0], nu2);

            GlobalPosition nu1(0);
            R.mv(connectVectors[1], nu1);

            // compute T, the area of quadrilateral made by normal vectors 'nu'
            GlobalPosition Rnu1(0);
            R.umv(nu1, Rnu1);
            Scalar T = fabs(nu2 * Rnu1);

            // Make new subVolume object
            SubVolume subVol;
            subVol.element = newElement;
            subVol.globalIdx = newElementGlobalIdx;
            subVol.localIsIndex[0] = IsIndicesInside[0];
            subVol.localIsIndex[1] = IsIndicesInside[1];
            subVol.contiPoints[0] = contiPoints[0];
            subVol.contiPoints[1] = contiPoints[1];
            subVol.nu[0] = nu1;
            subVol.nu[1] = nu2;
            subVol.T = T;

            // if we were on an internal flux or dirichlet boundary before and facet coupling is on,
            // adjust sub volume data -> i.e. contiPoint 1, nu 0 and T
            if ( (faceType == FaceTypes::InternalFluxFace || faceType == FaceTypes::InternalDirichletFace) && GET_PROP_VALUE(TypeTag, FacetCoupling) )
            {
                GlobalPosition newContiPoint2 = center;

                // get new connection vector from element to new conti point
                GlobalPosition newConnectVec = newContiPoint2;
                newConnectVec -= subVol.element.geometry().center();

                // calculate the new nu1
                GlobalPosition newNu1(0);
                R.mv(newConnectVec, newNu1);

                // compute new T, the area of quadrilateral made by normal vectors 'nu'
                GlobalPosition Rnu1(0);
                GlobalPosition oldNu2 = subVol.nu[1];
                R.umv(newNu1, Rnu1);
                Scalar newT = fabs(oldNu2 * Rnu1);

                // update data for the sub volume at hand
                subVol.contiPoints[1] = newContiPoint2;
                subVol.nu[0] = newNu1;
                subVol.T = newT;
            }

            // pass sub volume to interactionVolume container and increment counter
            interactionVolume->addSubVolume(subVol);
            interactionVolume->setElementToSubVolMap(newElementGlobalIdx,elementCounter);
            elementCounter++;

            // set the intersection iterator to the newly found intersection
            // and do the next face/volume
            nextIS = newIS;
        }

        // If volume has not been created yet we are on a boundary
        // of the counter clockwise rotation. Then we have to define the
        // face at hand as boundary and start the clockwise rotation using isIt14
        if (!volumeCreated)
        {
            // we definitely intersect the boundary here
            interactionVolume->setInteriorVolume(false);

            IntersectionGeometry isGeometry12 = nextIS->geometry();
            // get volume of face
            Scalar faceVol = isGeometry12.volume() / 2.0;
            // get outer normal vector of face
            GlobalPosition unitOuterNormal = nextIS->centerUnitOuterNormal();
            // get center of intersection
            GlobalPosition center12 = isGeometry12.center();

            BoundaryTypes bctypes;
            problem.boundaryTypesAtPos(bctypes, center12);

            // Make sub volume face object
            SubVolumeFace subFace;
            subFace.normal = unitOuterNormal;
            subFace.center = center12;
            subFace.area = faceVol;
            subFace.positiveSubVolIndex = elementCounter - 1;
            subFace.negativeSubVolIndex = -1;
            subFace.localIdxOnPositive = 0;
            subFace.localIdxOnNegative = -1;
            subFace.onBoundary = true;
            if (bctypes.isNeumann(phaseIdx))
            {
                subFace.faceType = FaceTypes::NeumannFace;
                interactionVolume->addInteriorFace(faceCounter);
            }
            else if (bctypes.isDirichlet(phaseIdx))
            {
                subFace.faceType = FaceTypes::DirichletFace;
                interactionVolume->addDirichletFace(faceCounter);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "MPFA-O method can only handle dirichlet/neumann boundaries so far! A different BC was tried to be applied.");
            subFace.boundaryTypes = bctypes;

            // pass sub volume face to interactionVolume container and increment counter
            interactionVolume->addSubVolumeFace(subFace);
            interactionVolume->setSubVolToFaceMap(subFace.positiveSubVolIndex, subFace.localIdxOnPositive, faceCounter);
            faceCounter++;

            // Now we start the clockwise rotation using isIt14
            IntersectionIterator nextIS2 = isIt14;
            int nextElementCounter = 0;
            while (!nextIS2->boundary())
            {
                IntersectionGeometry isGeometryNext = nextIS2->geometry();
                // get volume of face
                Scalar faceVolNext = isGeometryNext.volume() / 2.0;
                // get outer normal vector of face
                GlobalPosition unitOuterNormalNext = nextIS2->centerUnitOuterNormal();
                // get center of face
                GlobalPosition centerNext = isGeometryNext.center();

                // check if we are on an interior Boundary
                const int faceType = problem.isInternalBoundary(*nextIS2);

                // Make sub volume face object
                SubVolumeFace subFace;
                subFace.normal = unitOuterNormalNext;
                subFace.center = centerNext;
                subFace.area = faceVolNext;
                subFace.positiveSubVolIndex = nextElementCounter == 0 ? 0 : elementCounter - 1;
                subFace.negativeSubVolIndex = elementCounter;
                subFace.localIdxOnPositive = 1;
                subFace.localIdxOnNegative = 0;
                subFace.onBoundary = false;
                subFace.faceType = faceType;
                if (faceType == FaceTypes::InternalDirichletFace)
                    interactionVolume->addDirichletFace(faceCounter);
                else
                    interactionVolume->addInteriorFace(faceCounter);

                // pass sub volume face to interactionVolume container and increment counter
                interactionVolume->addSubVolumeFace(subFace);
                interactionVolume->setSubVolToFaceMap(subFace.positiveSubVolIndex, subFace.localIdxOnPositive, faceCounter);
                if (faceType != FaceTypes::InternalFluxFace && faceType != FaceTypes::InternalDirichletFace)
                    interactionVolume->setSubVolToFaceMap(subFace.negativeSubVolIndex, subFace.localIdxOnNegative, faceCounter);
                if ( (faceType == FaceTypes::InternalFluxFace || faceType == FaceTypes::InternalDirichletFace) && GET_PROP_VALUE(TypeTag, FacetCoupling) )
                {
                    GlobalPosition newContiPoint2 = subFace.center;
                    // get new connection vector from element to new conti point
                    GlobalPosition newConnectVec = newContiPoint2;
                    newConnectVec -= interactionVolume->getSubVolume(subFace.positiveSubVolIndex).element.geometry().center();

                    // calculate the new nu1
                    GlobalPosition newNu1(0);
                    R.mv(newConnectVec, newNu1);

                    // compute new T, the area of quadrilateral made by normal vectors 'nu'
                    GlobalPosition Rnu1(0);
                    GlobalPosition oldNu2 = interactionVolume->getSubVolume(subFace.positiveSubVolIndex).nu[1];
                    R.umv(newNu1, Rnu1);
                    Scalar newT = fabs(oldNu2 * Rnu1);

                    // update contiPoint for local index 0, nu for local index 1 and T
                    interactionVolume->updateSubVolumeData(subFace.positiveSubVolIndex, newContiPoint2, 1, newNu1, 0, newT);
                }
                faceCounter++;

                // if we are on an internal Flux boundary we have to duplicate the face swapping the sub volume indices
                if (faceType == FaceTypes::InternalFluxFace || faceType == FaceTypes::InternalDirichletFace)
                {
                    SubVolumeFace subFaceDouble;
                    subFaceDouble.normal = helpVector - unitOuterNormalNext;
                    subFaceDouble.center = centerNext;
                    subFaceDouble.area = faceVolNext;
                    subFaceDouble.positiveSubVolIndex = elementCounter;
                    subFaceDouble.negativeSubVolIndex = nextElementCounter == 0 ? 0 : elementCounter - 1;
                    subFaceDouble.localIdxOnPositive = 0;
                    subFaceDouble.localIdxOnNegative = 1;
                    subFaceDouble.onBoundary = false;
                    subFaceDouble.faceType = faceType;
                    if (faceType == FaceTypes::InternalDirichletFace)
                        interactionVolume->addDirichletFace(faceCounter);
                    else
                        interactionVolume->addInteriorFace(faceCounter);
                    // pass sub volume face to interactionVolume container and increment counter
                    interactionVolume->addSubVolumeFace(subFaceDouble);
                    // also set the mapping from the local face indices of the sub volumes to the face
                    interactionVolume->setSubVolToFaceMap(subFaceDouble.positiveSubVolIndex, subFaceDouble.localIdxOnPositive, faceCounter);
                    faceCounter++;
                }

                // obtain the information of the outside element
                Element newElement = nextIS2->outside();
                int newElementGlobalIdx = elementMapper.index(newElement);

                // find the intersection of the new element that contains the central node
                IntersectionIterator newIS;
                bool found = false;
                IntersectionIterator isNewEl = problem.gridView().ibegin(newElement);
                for(; isNewEl != problem.gridView().iend(newElement); ++isNewEl)
                {
                    if (found == true)
                        break;

                    for (int i = 0; i < isNewEl->geometry().corners(); ++i)
                    {
                        int localVertCorner = referenceElement.subEntity(isNewEl->indexInInside(), 1, i, dim);
                        int globalVertIdxcorner = vertexMapper.index( newElement.template subEntity<dim>(localVertCorner) );

                        if (globalVertIdxcorner == centralVertGlobalIdx && isNewEl->indexInInside() != nextIS2->indexInOutside())
                        {
                            newIS = isNewEl;  found = true;
                            break;
                        }
                    }
                }
                if (found == false)
                    DUNE_THROW(Dune::InvalidStateException, "Couldn't find facet of the neighbouring element that contains the same node!");

                // copy the two intersection indices in a container and calculate
                // continuity points. Boolean swap used as an error indicator
                std::vector<int> IsIndicesInside(2);
                IsIndicesInside[0] = nextIS2->indexInOutside();
                IsIndicesInside[1] = newIS->indexInInside();

                std::vector<GlobalPosition> connectVectors(2);
                std::vector<GlobalPosition> contiPoints(2);
                bool swap = calcConnectionVectors(newElement, IsIndicesInside, centralVertPosition, connectVectors, contiPoints);
                if (swap)
                    DUNE_THROW(Dune::InvalidStateException, "This shouldn't happen in mpfao2dmanager.hh, l. 448");

                // compute normal vectors nu1,nu2
                GlobalPosition nu2(0);
                R.mv(helpVector - connectVectors[0], nu2);

                GlobalPosition nu1(0);
                R.mv(connectVectors[1], nu1);

                // compute T, the area of quadrilateral made by normal vectors 'nu'
                GlobalPosition Rnu1(0);
                R.umv(nu1, Rnu1);
                Scalar T = fabs(nu2 * Rnu1);

                // Make new subVolume object
                SubVolume subVol;
                subVol.element = newElement;
                subVol.globalIdx = newElementGlobalIdx;
                subVol.localIsIndex[0] = IsIndicesInside[0];
                subVol.localIsIndex[1] = IsIndicesInside[1];
                subVol.contiPoints[0] = contiPoints[0];
                subVol.contiPoints[1] = contiPoints[1];
                subVol.nu[0] = nu1;
                subVol.nu[1] = nu2;
                subVol.T = T;

                // if we were on an internal dirichlet face before and facet coupling is on,
                // adjust sub volume data -> i.e. contiPoint 1, nu 0 and T
                if ( (faceType == FaceTypes::InternalFluxFace || faceType == FaceTypes::InternalDirichletFace) && GET_PROP_VALUE(TypeTag, FacetCoupling) )
                {
                    GlobalPosition newContiPoint1 = centerNext;

                    // get new connection vector from element to new conti point
                    GlobalPosition newConnectVec = newContiPoint1;
                    newConnectVec -= subVol.element.geometry().center();

                    // calculate the new nu2
                    GlobalPosition newNu2(0);
                    R.mv(helpVector - newConnectVec, newNu2);

                    // compute new T, the area of quadrilateral made by normal vectors 'nu'
                    GlobalPosition Rnu1(0);
                    GlobalPosition oldNu1 = subVol.nu[0];
                    R.umv(oldNu1, Rnu1);
                    Scalar newT = fabs(newNu2 * Rnu1);

                    // update data for the sub volume at hand
                    subVol.contiPoints[0] = newContiPoint1;
                    subVol.nu[1] = newNu2;
                    subVol.T = newT;
                }

                // pass sub volume to interactionVolume container and increment counter
                interactionVolume->addSubVolume(subVol);
                interactionVolume->setElementToSubVolMap(newElementGlobalIdx,elementCounter);
                elementCounter++; nextElementCounter++;

                // set our intersection iterator to the newly found intersection
                // and do the next face/volume
                nextIS2 = newIS;
            }

            // If we are here, we need to define the boundary face
            // of the counterclockwise rotation
            IntersectionGeometry isGeometry2 = nextIS2->geometry();
            Scalar faceVol2 = isGeometry2.volume() / 2.0;
            // get outer normal vector of face
            GlobalPosition unitOuterNormal2 = nextIS2->centerUnitOuterNormal();
            // get center of the intersection
            GlobalPosition center2 = isGeometry2.center();

            BoundaryTypes bctypes2;
            problem.boundaryTypesAtPos(bctypes2, center2);

            // Make sub volume face object
            SubVolumeFace subFace2;
            subFace2.normal = unitOuterNormal2;
            subFace2.center = center2;
            subFace2.area = faceVol2;
            subFace2.positiveSubVolIndex = nextElementCounter == 0 ? 0 : elementCounter - 1;
            subFace2.negativeSubVolIndex = -1;
            subFace2.localIdxOnPositive = 1;
            subFace2.localIdxOnNegative = -1;
            subFace2.onBoundary = true;
            if (bctypes2.isNeumann(phaseIdx))
            {
                subFace2.faceType = FaceTypes::NeumannFace;
                interactionVolume->addInteriorFace(faceCounter);
            }
            else if (bctypes2.isDirichlet(phaseIdx))
            {
                subFace2.faceType = FaceTypes::DirichletFace;
                interactionVolume->addDirichletFace(faceCounter);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "MPFA-O method can only handle dirichlet/neumann boundaries so far! A different BC was tried to be applied.");
            subFace2.boundaryTypes = bctypes2;

            // pass sub volume face to interactionVolume container and increment counter
            interactionVolume->addSubVolumeFace(subFace2);
            interactionVolume->setSubVolToFaceMap(subFace2.positiveSubVolIndex, subFace2.localIdxOnPositive, faceCounter);
            faceCounter++;
        }

        // pass the number of sub volumes and subfaces found in the interaction volume
        interactionVolume->setNumberOfSubVols(elementCounter);
        interactionVolume->setNumberOfSubFaces(faceCounter);

        // set up transmissivity matrix and store in interaction volume
        interactionVolume->calculateTransmissivityMatrix(problem);
        // uncomment lines below to have the interactionvolume and its transmissivity matrix printed out
        //interactionVolume->printInteractionVolume();
        //interactionVolume->printTransmissivityMatrix();
    } // end of function


    // fills a storage vector with the element pointers of the stencil
    // assumes that the first element of the input container
    // contains a pointer to the element of which the stencil is desired
    static int getElementsOfStencil(std::vector<Element> &container, const Problem &problem)
    {
        int numNeighbors = 1;

        // Loop over the nodes of the element and the corresponding interaction volumes
        // to fill the the container with info about the elements of the stencil
        for (int nodeIdx = 0; nodeIdx < container[0].subEntities(dim); nodeIdx++)
        {
            std::vector<Element> tmpElements;

            int boundaryLayer = problem.model().interactionVolumeContainer().isBoundaryVolume(container[0], /*dummy*/0, nodeIdx);

            if (boundaryLayer == BoundaryLayers::interior)
                problem.model().getInsideInteractionVolume(container[0], /*dummy*/0, nodeIdx)->passElementsInRegion(tmpElements);
            else
                problem.model().getBoundaryInteractionVolume(container[0], /*dummy*/0, nodeIdx, /*dummy*/0)->passElementsInRegion(tmpElements);

            // insert new elements in neighbor container
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
        const GlobalPosition isCenter1 = element.template subEntity<1>(IsIndicesInsideElement[0]).geometry().center();
        const GlobalPosition isCenter2 = element.template subEntity<1>(IsIndicesInsideElement[1]).geometry().center();

        GlobalPosition connectVec1 = corner;
        connectVec1 -= isCenter1;
        connectVec1 *= continuityPoint;
        GlobalPosition connectVec2 = corner;
        connectVec2 -= isCenter2;
        connectVec2 *= continuityPoint;

        // calculate the positions on which flux and pressure continuity is demanded
        GlobalPosition contiPoint1 = isCenter1;
        contiPoint1 += connectVec1;
        GlobalPosition contiPoint2 = isCenter2;
        contiPoint2 += connectVec2;

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
            DUNE_THROW(Dune::InvalidStateException, "Error! Dot Product = 0 in \"makeRightHandedSystem \" in mpfaointeractionvolume2dfiller.hh");

        return swap;
    }
};
}
#endif
