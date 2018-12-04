// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
#ifndef DUMUX_ONEP_TUBES_TEST_PROBLEM_HH
#define DUMUX_ONEP_TUBES_TEST_PROBLEM_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/reorderingdofmapper.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>


#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class TubesTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct TubesTest { using InheritsFrom = std::tuple<OneP>; };
struct TubesTestCCTpfa { using InheritsFrom = std::tuple<TubesTest, CCTpfaModel>; };
struct TubesTestBox { using InheritsFrom = std::tuple<TubesTest, BoxModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_FOAMGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TubesTest> { using type = Dune::FoamGrid<1, 3>; };
#endif

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (cc)
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::TubesTestCCTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using ElementMapper = ReorderingDofMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (box)
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::TubesTestBox>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = ReorderingDofMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TubesTest> { using type = TubesTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TubesTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TubesTestSpatialParams<FVGridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TubesTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};
} // end namespace Properties

/*!
 * Miminmimal setting for root systems,
 * constant soil, constant transpiration from input file
 */
template <class TypeTag>
class TubesTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    // Grid and world dimension
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box };
    enum {
        bcDirichlet = 0, bcNeumann = 1
    };

public:
    TubesTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
:
        ParentType(fvGridGeometry)
    {
        soilP_ = toPa_(getParam<Scalar>("RootSystem.Soil.P"));
        try {
            collarP_ = toPa_(getParam<Scalar>("RootSystem.Collar.P"));
            type_ = bcDirichlet;
        } catch (std::exception &e) {
            trans_ = toPa_(getParam<Scalar>("RootSystem.Collar.Transpirationb")) / 100. * 1000. / (3600 * 24); // cm/day *rho -> [ kg / (m^2 \cdot s)]
            type_ = bcNeumann;
        }
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    {
        return 273.15 + 37.0;  // Body temperature
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
        const SubControlVolume &scv, const ElementSolution& elemSol) const
    {
        const auto a = this->spatialParams().radius(element);
        // std::cout << "extrusion factor radius " << a;
        return M_PI * a * a;
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // default
        if (onUpperBoundary_(pos)) { // root collar
            if (type_ == bcDirichlet) {
                bcTypes.setAllDirichlet();
            } else {
                bcTypes.setAllNeumann();
            }
        } else { // for all other (i.e. root tips)
            std::cout << "Set noflux bc \n";
            bcTypes.setAllNeumann();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const
    {
        return PrimaryVariables(collarP_);
    }

    /*
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumannAtPos(const GlobalPosition &pos) const {
        if (onUpperBoundary_(pos)) {
            return NumEqVector(trans_);
        } else {
            return NumEqVector(0.);
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    NumEqVector source(const Element &element,
        const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const
    {
        NumEqVector values;
        auto params = this->spatialParams();
        Scalar a = params.radius(element); // root radius (m)
        Scalar kr = params.radialConductivity(element); //  radial conductivity (m^2 s / kg)
        Scalar phx = elemVolVars[0].pressure(); // Pa, kg/m/s^2
        Scalar phs = soilP_; // Pa, kg/m/s^2
        // std::cout << "problem " << phs << ", " << a << ", " << kr * 1.e9 << "\n";
        values[conti0EqIdx] = kr * 2 * a * M_PI * (phs - phx); // m^3/s
        values[conti0EqIdx] /= (a * a * M_PI); // 1/s
        values[conti0EqIdx] *= rho_; // (kg/s/m^3)
        std::cout << values << "\n";
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(soilP_);
    }

    // \}

private:

    Scalar soilP_;
    int type_ = 0;
    Scalar collarP_ = 0.; // pa
    Scalar trans_ = 0.; // m/s

    Scalar toPa_(Scalar ph) const { // cm -> Pa
        return pRef_ + ph / 100. * rho_ * g_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const { // on root collar
        return globalPos[dimWorld - 1] > this->fvGridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa
    static constexpr Scalar eps_ = 1e-8;
};

} //end namespace Dumux

#endif
