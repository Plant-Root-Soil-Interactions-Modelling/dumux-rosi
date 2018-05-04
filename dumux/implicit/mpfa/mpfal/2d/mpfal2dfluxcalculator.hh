#ifndef DUMUX_MPFAL2DFLUXCALCULATOR_HH
#define DUMUX_MPFAL2DFLUXCALCULATOR_HH

#include <dumux/common/math.hh>
#include <eigen3/Eigen/Dense>

namespace Dumux{

template<class TypeTag>
class MpfaL2DFluxCalculator
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolume) BoundaryInteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename InteractionVolume::TransmissivityMatrix TransMatrix;
    typedef typename BoundaryInteractionVolume::DynamicVector FluxVector;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:

    static Scalar calculateScvFaceFlux(const Problem &problem, const Element &element, int localIsIndex,
                                const InteractionVolume &interactionVolume, const ElementVolumeVariables &elemVolVars,
                                const FVElementGeometry &fvGeometry, int phaseIdx, int &regionFaceIdx, FluxVector &neumannFluxes)
    {
        const TransMatrix &T = interactionVolume.getTransMatrix();

        // Set up vector with the primary variable values
        Dune::FieldVector<Scalar, 3> piezoHeads(0);
        Dune::FieldVector<Scalar, 2> fluxes(0);

        // find the face index of the facet in the interaction region
        // and determine if the flux has to be multiplied by -1
        bool switchSign = false;
        if (interactionVolume.secondTriangle() && interactionVolume.getSubVolIdx(element) != 0)
            switchSign = true;
        else if (!interactionVolume.secondTriangle() && interactionVolume.getSubVolIdx(element) != 0)
            switchSign = true;
        int faceIndexInRegion = interactionVolume.getFluxFaceIdx();
        regionFaceIdx = faceIndexInRegion;

        // set up piezometric head vector entries of the sub volumes (3 sub volumes in 2D)
        for (int head = 0; head < 3; head++)
        {
            const Element &subVol = interactionVolume.getSubVolume(head);
            // get global index of sub volume element
            int globalElemIdx = problem.elementMapper().index(subVol);

            int idxInStencil = fvGeometry.findElementInStencil(globalElemIdx);
            Scalar pressure = elemVolVars[idxInStencil].fluidState().pressure(phaseIdx);

            // turn pressure into piezometric head
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // ask for the gravitational acceleration at the given point
                const GlobalPosition subVolCenter = subVol.geometry().center();
                GlobalPosition g( problem.gravityAtPos(subVolCenter) );
                const Scalar density = elemVolVars[idxInStencil].fluidState().density(phaseIdx);
                // make gravity acceleration a force
                Scalar f = g * subVolCenter;
                f *= density;
                // calculate the final potential
                pressure -= f;
            }
            // set value in the piezometric heads vector
            piezoHeads[head] = pressure;
        }

        T.umv(piezoHeads, fluxes);

        if (switchSign)
            return -fluxes[faceIndexInRegion];
        else
            return fluxes[faceIndexInRegion];
    }
};
}
#endif
