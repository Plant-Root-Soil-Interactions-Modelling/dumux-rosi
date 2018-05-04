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
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Intersection computation relies on the dune grid-glue backend.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_CIRCLESOURCES_HH
#define DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_CIRCLESOURCES_HH

#include <dumux/multidimension/embeddedcoupling/cellcentered/bboxtreecouplingmanager.hh>

namespace Dumux
{

/*!
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 * The pressure for the source term is computed as an integral over the circle describing
 * the surface of the one-dimensional tube. This exact determination of the bulk pressure
 * is necessary for convergence studies.
 */
template<class TypeTag>
class CCBBoxTreeEmbeddedCouplingManagerCircleSources : public CCBBoxTreeEmbeddedCouplingManager<TypeTag>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ParentType = Dumux::CCBBoxTreeEmbeddedCouplingManager<TypeTag>;

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkGrid = typename GET_PROP_TYPE(BulkProblemTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid);

    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    using PointSourceData = Dumux::PointSourceDataCircleAverage<TypeTag>;

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    /*!
     * \brief Constructor
     */
    CCBBoxTreeEmbeddedCouplingManagerCircleSources(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : ParentType(bulkProblem, lowDimProblem) {}

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element)
     * \param order The order of the quadrature rule for integration of sources over an intersection
     * \param verbose If the point source computation is verbose
     */
    void computePointSourceData(unsigned int order = 1, bool verbose = false)
    {
        // Initialize the bulk bounding box tree
        auto bulkTree = this->bulkProblem().boundingBoxTree();

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // iterate over all lowdim elements
        for (const auto& lowDimElement : elements(this->lowDimGridView()))
        {
            // get the Gaussian quadrature rule for the low dim element
            auto lowDimGeometry = lowDimElement.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(lowDimGeometry.type(), order);

            const unsigned int lowDimElementIdx = this->lowDimProblem().elementMapper().index(lowDimElement);

            // apply the Gaussian quadrature rule and define point sources at each quadrature point
            // note that the approximation is not optimal if
            // (a) the one-dimensional elements are too large,
            // (b) whenever a one-dimensional element is split between two or more elements,
            // (c) when gradients of important quantities in the three-dimensional domain are large.

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // global position of the quadrature point
                const auto globalPos = lowDimGeometry.global(qp.position());

                auto bulkElementIndices = bulkTree.computeEntityCollisions(globalPos);

                // do not add a point source if the qp is outside of the 3d grid
                // this is equivalent to having a source of zero for that qp
                if (bulkElementIndices.empty())
                    continue;

                //////////////////////////////////////////////////////////
                // get circle average connectivity and interpolation data
                //////////////////////////////////////////////////////////

                auto numIp = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, MultiDimension, NumCircleSegments);
                const auto radius = this->lowDimProblem().spatialParams().radius(lowDimElementIdx);
                const auto normal = lowDimGeometry.corner(1)-lowDimGeometry.corner(0);
                auto length = 2*M_PI*radius/numIp;

                auto circlePoints = this->getCirclePoints_(globalPos, normal, radius, numIp);
                std::vector<Scalar> circleIpWeight;
                std::vector<unsigned int> circleStencil;
                for (const auto& p : circlePoints)
                {
                    auto bulkElementIndices = bulkTree.computeEntityCollisions(p);
                    if (bulkElementIndices.empty())
                        continue;

                    for (auto bulkElementIdx : bulkElementIndices)
                    {
                        const auto id = idCounter_++;
                        const auto ie = lowDimGeometry.integrationElement(qp.position())*length;
                        const auto qpweight = 1;

                        this->bulkPointSources().emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({bulkElementIdx}));
                        this->bulkPointSources().back().setEmbeddings(bulkElementIndices.size());

                        // compute the coupling stencils
                        this->bulkCouplingStencils()[bulkElementIdx].push_back(lowDimElementIdx);

                        PointSourceData psData;
                        psData.addLowDimInterpolation(lowDimElementIdx);
                        psData.addBulkInterpolation(bulkElementIdx);

                        // publish point source data in the global vector
                        pointSourceData_.push_back(psData);

                        circleStencil.push_back(bulkElementIdx);
                        circleIpWeight.push_back(length);
                    }
                }

                // export low dim circle stencil
                lowDimCircleStencils_[lowDimElementIdx].insert(lowDimCircleStencils_[lowDimElementIdx].end(),
                                                               circleStencil.begin(),
                                                               circleStencil.end());

                for (auto bulkElementIdx : bulkElementIndices)
                {
                    const auto id = idCounter_++;
                    const auto ie = lowDimGeometry.integrationElement(qp.position());
                    const auto qpweight = qp.weight();

                    this->lowDimPointSources().emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({lowDimElementIdx}));
                    this->lowDimPointSources().back().setEmbeddings(bulkElementIndices.size());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);
                    // add data needed to compute integral over the circle
                    psData.addCircleInterpolation(circleIpWeight, circleStencil);

                    // publish point source data in the global vector
                    pointSourceData_.push_back(psData);
                }
            }
        }

        // make the circle stencils unique
        for (auto&& stencil : bulkCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        for (auto&& stencil : lowDimCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // coupling stencils
        for (auto&& stencil : this->bulkCouplingStencils())
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // the low dim coupling stencil
        for (auto&& stencil : this->lowDimCouplingStencils())
        {
            // add the circle stencil to the coupling stencil
            stencil.second.insert(stencil.second.end(),
                                  lowDimCircleStencils_.at(stencil.first).begin(),
                                  lowDimCircleStencils_.at(stencil.first).end());
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // compute the bulk stencil (like in the fv element geometry)
        for (const auto& element : elements(this->bulkGridView()))
        {
            // the element itself is first
            const unsigned int eIdx = this->bulkProblem().elementMapper().index(element);
            this->bulkStencils().insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(this->bulkGridView(), element))
                if (intersection.neighbor())
                    this->bulkStencils().at(eIdx).push_back(this->bulkProblem().elementMapper().index(intersection.outside()));
        }

        // ...and afterwards add the bulk circle stencil so that the order doesn't get messed up.
        for (auto&& stencil : bulkCircleStencils_)
        {
            auto&& bulkStencil = this->bulkStencils().at(stencil.first);
            // only insert elements that are not already in the stencil
            for (auto&& eIdx : stencil.second)
                if (std::find(bulkStencil.begin(), bulkStencil.end(), eIdx) == bulkStencil.end())
                    bulkStencil.push_back(eIdx);
        }

        // compute the low dimensional stencil (like in the fv geometry)
        for (const auto& element : elements(this->lowDimGridView()))
        {
            // the element itself is first
            const unsigned int eIdx = this->lowDimProblem().elementMapper().index(element);
            this->lowDimStencils().insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(this->lowDimGridView(), element))
                if (intersection.neighbor())
                    this->lowDimStencils().at(eIdx).push_back(this->lowDimProblem().elementMapper().index(intersection.outside()));
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(unsigned int id) const
    {
        return pointSourceData_[id];
    }

    //! Return data for a bulk point source with the identifier id
    const BulkPrimaryVariables& bulkPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateBulk(this->bulkProblem().model().curSol());
        return data.bulkPriVars();
    }

    //! Return data for a low dim point source with the identifier id
    const LowDimPrimaryVariables& lowDimPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateLowDim(this->lowDimProblem().model().curSol());
        return data.lowDimPriVars();
    }

    //! Return a reference to the bulk problem
    Scalar radius(unsigned int id) const
    {
        const auto& data = pointSourceData_[id];
        return this->lowDimProblem().spatialParams().radius(data.lowDimElementIdx());
    }
    // \}

    //! Clear all internal data members
    void clear()
    {
        ParentType::clear();
        pointSourceData_.clear();
        bulkCircleStencils_.clear();
        lowDimCircleStencils_.clear();
        idCounter_ = 0;
    }

private:
    std::vector<GlobalPosition> getCirclePoints_(const GlobalPosition& center,
                                                 const GlobalPosition& normal,
                                                 const Scalar radius,
                                                 const unsigned int numIp = 20)
    {
        const Scalar eps = 1.5e-7;
        static_assert(dimWorld == 3, "Only implemented for world dimension 3");

        std::vector<GlobalPosition> points(numIp);

        // make sure n is a unit vector
        auto n = normal;
        n /= n.two_norm();

        // caculate a vector u perpendicular to n
        GlobalPosition u;
        if (std::abs(n[0]) < eps && std::abs(n[1]) < eps)
            if (std::abs(n[2]) < eps)
                DUNE_THROW(Dune::MathError, "The normal vector has to be non-zero!");
            else
                u = {0, 1, 0};
        else
            u = {-n[1], n[0], 0};

        u *= radius/u.two_norm();

        // the circle parameterization is p(t) = r*cos(t)*u + r*sin(t)*(n x u) + c
        auto tangent = crossProduct(u, n);
        tangent *= radius/tangent.two_norm();

        // the parameter with an offset
        Scalar t = 0 + 0.1;
        // insert the vertices
        for (unsigned int i = 0; i < numIp; ++i)
        {
            points[i] = GlobalPosition({u[0]*std::cos(t) + tangent[0]*std::sin(t) + center[0],
                                        u[1]*std::cos(t) + tangent[1]*std::sin(t) + center[1],
                                        u[2]*std::cos(t) + tangent[2]*std::sin(t) + center[2]});

            t += 2*M_PI/numIp;

            // periodic t
            if(t > 2*M_PI)
                t -= 2*M_PI;
        }

        return std::move(points);
    }

    // id generator for point sources
    unsigned int idCounter_;

    // special point source data for the circle average
    mutable std::vector<PointSourceData> pointSourceData_;

    // circle stencils
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCircleStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCircleStencils_;
};

} // namespace Dumux

#endif
