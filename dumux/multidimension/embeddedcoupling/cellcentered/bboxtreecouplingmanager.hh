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
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH
#define DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/multidimension/glue/glue.hh>
#include <dumux/multidimension/embeddedcoupling/pointsourcedata.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

/*!
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class TypeTag>
class CCBBoxTreeEmbeddedCouplingManager
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkPointSource = typename GET_PROP_TYPE(BulkProblemTypeTag, PointSource);
    using LowDimPointSource = typename GET_PROP_TYPE(LowDimProblemTypeTag, PointSource);

    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);

    using BulkVolumeVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, VolumeVariables);
    using LowDimVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, VolumeVariables);

    using BulkFVElementGeometry = typename GET_PROP_TYPE(BulkProblemTypeTag, FVElementGeometry);
    using LowDimFVElementGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVElementGeometry);

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

    using PointSourceData = Dumux::PointSourceData<TypeTag>;

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

public:

    /*!
     * \brief Constructor
     */
    CCBBoxTreeEmbeddedCouplingManager(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblem_(bulkProblem), lowDimProblem_(lowDimProblem)
    {
        // Check if we are using the cellcentered method in both domains
        static_assert(!bulkIsBox && !lowDimIsBox, "Using the cell-centered coupling manager for problems using box discretization!");
        static_assert(lowDimDim == 1, "The bounding box coupling manager only works with one-dimensional low-dim grids");
    }

    /*!
     * \brief Methods to be accessed by the coupled problem
     */
    // \{

    /*!
     * \brief Called by the coupled problem before
     *        initializing the subproblems / models
     */
    void preInit()
    {
        asImp_().computePointSourceData();
    }

    /*!
     * \brief Called by the coupled problem after
     *        initializing the subproblems / models
     */
    void postInit()
    {}

    /*!
     * \brief The bulk coupling stencil, i.e. which low-dimensional dofs
     *        the given bulk element dof depends on.
     */
    const std::vector<unsigned int>& couplingStencil(const BulkElement& element) const
    {
        const unsigned int eIdx = bulkProblem().elementMapper().index(element);
        if (bulkCouplingStencils_.count(eIdx))
            return bulkCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The low dim coupling stencil, i.e. which bulk dofs
     *        the given low dimensional element dof depends on.
     */
    const std::vector<unsigned int>& couplingStencil(const LowDimElement& element) const
    {
        const unsigned int eIdx = lowDimProblem().elementMapper().index(element);
        if (lowDimCouplingStencils_.count(eIdx))
            return lowDimCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The bulk stencil, i.e. which bulk dofs
     *        the given bulk element dof depends on.
     */
    const std::vector<unsigned int>& stencil(const BulkElement& element) const
    {
        const unsigned int eIdx = bulkProblem().elementMapper().index(element);
        assert(bulkStencils_.count(eIdx)); // should never be empty
        return bulkStencils_.at(eIdx);
    }

    /*!
     * \brief The low dim stencil, i.e. which low-dimensional dofs
     *        the given low-dimensional element dof depends on.
     */
    const std::vector<unsigned int>& stencil(const LowDimElement& element) const
    {
        const unsigned int eIdx = lowDimProblem().elementMapper().index(element);
        assert(lowDimStencils_.count(eIdx)); // should never be empty
        return lowDimStencils_.at(eIdx);
    }

    // \}

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
        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // a vector for the volumes the low-dim domain occupies in the bulk domain if it were full-dimensional
        lowDimVolumeInBulkElement_.resize(this->bulkGridView().size(0));
        // a vector for the total low-dim tube surface in the bulk domain
        lowDimSurfaceInBulkElement_.resize(this->bulkGridView().size(0));

        // intersect the bounding box trees
        Dumux::CCMultiDimensionGlue<TypeTag> glue(bulkProblem(), lowDimProblem());
        glue.build();

        for (const auto& is : intersections(glue))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            using Quad = Dune::QuadratureRule<Scalar, lowDimDim>;
            const Quad &quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);

            const unsigned int lowDimElementIdx = lowDimProblem().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = lowDimProblem().spatialParams().radius(lowDimElementIdx);
            for (unsigned int outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);
                if (lowDimDim == 1)
                {
                    lowDimVolumeInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*M_PI*radius*radius;
                    lowDimSurfaceInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*2*M_PI*radius;
                    lowDimStencilLength_[std::make_pair(lowDimElementIdx, bulkElementIdx)] = intersectionGeometry.volume();
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Volume computation for low dimensional domain of dim " << lowDimDim);
            }

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (unsigned int outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
                {
                    const auto& outside = is.outside(outsideIdx);
                    const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    bulkPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({bulkElementIdx}));
                    bulkPointSources_.back().setEmbeddings(is.neighbor(0));
                    lowDimPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({lowDimElementIdx}));
                    lowDimPointSources_.back().setEmbeddings(is.neighbor(0));

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);

                    // publish point source data in the global vector
                    pointSourceData_.push_back(psData);

                    // compute the coupling stencils
                    bulkCouplingStencils_[bulkElementIdx].push_back(lowDimElementIdx);
                    // add this bulk element to the low dim coupling stencil
                    lowDimCouplingStencils_[lowDimElementIdx].push_back(bulkElementIdx);
                }
            }
        }

        // coupling stencils
        for (auto&& stencil : bulkCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }
        // the low dim coupling stencil
        for (auto&& stencil : lowDimCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // compute the bulk stencil (like in the fv element geometry)
        for (const auto& element : elements(bulkGridView()))
        {
            // the element itself is first
            const unsigned int eIdx =  bulkProblem().elementMapper().index(element);
            bulkStencils_.insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(bulkGridView(), element))
                if (intersection.neighbor())
                    bulkStencils_.at(eIdx).push_back(bulkProblem().elementMapper().index(intersection.outside()));
        }

        // compute the low dimensional stencil (like in the fv geometry)
        for (const auto& element : elements(lowDimGridView()))
        {
            // the element itself is first
            const unsigned int eIdx =  lowDimProblem().elementMapper().index(element);
            lowDimStencils_.insert({eIdx, {eIdx}});

            for (const auto& intersection : intersections(lowDimGridView(), element))
                if (intersection.neighbor())
                    lowDimStencils_.at(eIdx).push_back(lowDimProblem().elementMapper().index(intersection.outside()));
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the bulk problem
    Scalar radius(unsigned int id) const
    {
        const auto& data = pointSourceData_[id];
        return lowDimProblem().spatialParams().radius(data.lowDimElementIdx());
    }

    //! The volume the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolume(const BulkElement& element) const
    {
        auto eIdx = bulkProblem().elementMapper().index(element);
        return lowDimVolumeInBulkElement_[eIdx];
    }

    //! The total low-dim tube surface in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimSurface(const BulkElement& element) const
    {
        auto eIdx = bulkProblem().elementMapper().index(element);
        return lowDimSurfaceInBulkElement_[eIdx];
    }

    //! The total low-dim tube surface in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimSurface(const int eIdx) const
    {
        return lowDimSurfaceInBulkElement_[eIdx];
    }

    //! The volume fraction the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolumeFraction(const BulkElement& element) const
    {
        auto totalVolume = element.geometry().volume();
        return asImp_().lowDimVolume(element) / totalVolume;
    }

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(unsigned int id) const
    {
        return pointSourceData_[id];
    }

    //! Return a reference to the bulk problem
    const BulkProblem& bulkProblem() const
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    const LowDimProblem& lowDimProblem() const
    {
        return lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    BulkProblem& bulkProblem()
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    LowDimProblem& lowDimProblem()
    {
        return lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    const BulkGridView& bulkGridView() const
    {
        return bulkProblem().gridView();
    }

    //! Return a reference to the low dimensional problem
    const LowDimGridView& lowDimGridView() const
    {
        return lowDimProblem().gridView();
    }

    //! Return a reference to the point sources
    const std::vector<BulkPointSource>& bulkPointSources() const
    {
        return bulkPointSources_;
    }

    //! Return a reference to the low dimensional problem
    const std::vector<LowDimPointSource>& lowDimPointSources() const
    {
        return lowDimPointSources_;
    }

    //! Return data for a bulk point source with the identifier id
    const BulkPrimaryVariables& bulkPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateBulk(bulkProblem().model().curSol());
        return data.bulkPriVars();
    }

    //! Return data for a low dim point source with the identifier id
    const LowDimPrimaryVariables& lowDimPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateLowDim(lowDimProblem().model().curSol());
        return data.lowDimPriVars();
    }

    //! Compute bulk volume variables for the identifier id
    BulkVolumeVariables bulkVolVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateBulk(bulkProblem().model().curSol());
        const auto& priVars = data.bulkPriVars();

        auto element = bulkProblem().boundingBoxTree().entity(data.bulkElementIdx());
        BulkFVElementGeometry fvGeometry;
        fvGeometry.update(bulkGridView(), element);

        BulkVolumeVariables volVars;
        volVars.update(priVars,
                       bulkProblem(),
                       element,
                       fvGeometry,
                       0, false);

        return volVars;
    }

    //! Compute lowDim volume variables for the identifier id
    LowDimVolumeVariables lowDimVolVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        data.interpolateLowDim(lowDimProblem().model().curSol());
        const auto& priVars = data.lowDimPriVars();

        auto element = lowDimProblem().boundingBoxTree().entity(data.lowDimElementIdx());
        LowDimFVElementGeometry fvGeometry;
        fvGeometry.update(lowDimGridView(), element);

        LowDimVolumeVariables volVars;
        volVars.update(priVars,
                       lowDimProblem(),
                       element,
                       fvGeometry,
                       0, false);

        return volVars;
    }
    // \}

    //! Clear all internal data members
    void clear()
    {
        bulkPointSources_.clear();
        lowDimPointSources_.clear();
        pointSourceData_.clear();
        bulkCouplingStencils_.clear();
        lowDimCouplingStencils_.clear();
        bulkStencils_.clear();
        lowDimStencils_.clear();
        idCounter_ = 0;
    }
    //! Return reference to low dim coupling stencil member
    const std::unordered_map<unsigned int, std::vector<unsigned int> >& lowDimCouplingStencils() const
    { return lowDimCouplingStencils_; }

    //! Initialize Coupling sources: set value to 0 and reuseCouplingSource to false
    void initializeCouplingSources()
    {   for (auto&& stencil : lowDimCouplingStencils())
        {
            for (auto&& second : stencil.second)
            {
                int rootIdx = stencil.first;
                int soilIdx = second;
                couplingSources_[std::make_pair(rootIdx, soilIdx)] = 0.;
                reuseCouplingSources_[std::make_pair(rootIdx, soilIdx)] = false;
            }
        }
    }

    //! Set coupling sources to values
    void setCouplingSources(std::pair<unsigned int, unsigned int> rootIdx_soilIdx, BulkPrimaryVariables value ) const
    {
        couplingSources_[rootIdx_soilIdx] = value;
    }

    //! Return values of coupling sources
    BulkPrimaryVariables couplingSources(std::pair<unsigned int, unsigned int> rootIdx_soilIdx) const
    {
        return couplingSources_[rootIdx_soilIdx];
    }

    //! Set True to reuse the values of coupling sources
    void setReuseCouplingSources(std::pair<unsigned int, unsigned int> rootIdx_soilIdx) const
    {
        reuseCouplingSources_[rootIdx_soilIdx] = true;
    }

    //! Return bool to decide whether to use the coupling sources
    bool reuseCouplingSources(std::pair<unsigned int, unsigned int> rootIdx_soilIdx) const
    {
        return reuseCouplingSources_[rootIdx_soilIdx];
    }

    //! Set reuseCouplingSources to false
    void resetCouplingSources() const
    {
        for (auto&& stencil : lowDimCouplingStencils())
        {
            for (auto&& second : stencil.second)
            {
                int rootIdx = stencil.first;
                int soilIdx = second;
                reuseCouplingSources_[std::make_pair(rootIdx, soilIdx)] = false;
            }
        }
    }

    //! Return total values of coupling sources
    BulkPrimaryVariables totalCouplingSources() const
    {
        BulkPrimaryVariables totalCouplingSources_;
        for (auto&& stencil : lowDimCouplingStencils())
        {
            for (auto&& second : stencil.second)
            {
                int rootIdx = stencil.first;
                int soilIdx = second;
                totalCouplingSources_ += couplingSources_[std::make_pair(rootIdx, soilIdx)];
            }
        }
        return totalCouplingSources_;
    }

    //! Return the length of intersected low-dim segment
    Scalar lowDimStencilLength(std::pair<unsigned int, unsigned int> rootIdx_soilIdx) const
    {
        return lowDimStencilLength_[rootIdx_soilIdx];
    }

protected:
    //! Return reference to point source data vector member
    std::vector<PointSourceData>& pointSourceData()
    { return pointSourceData_; }

    //! Return reference to bulk point sources
    std::vector<BulkPointSource>& bulkPointSources()
    { return bulkPointSources_; }

    //! Return reference to low dim point sources
    std::vector<LowDimPointSource>& lowDimPointSources()
    { return lowDimPointSources_; }

    //! Return reference to bulk coupling stencil member
    std::unordered_map<unsigned int, std::vector<unsigned int> >& bulkCouplingStencils()
    { return bulkCouplingStencils_; }

    //! Return reference to low dim coupling stencil member
    std::unordered_map<unsigned int, std::vector<unsigned int> >& lowDimCouplingStencils()
    { return lowDimCouplingStencils_; }

    //! Return reference to bulk stencil member
    std::unordered_map<unsigned int, std::vector<unsigned int> >& bulkStencils()
    { return bulkStencils_; }

    //! Return reference to low dim stencil member
    std::unordered_map<unsigned int, std::vector<unsigned int> >& lowDimStencils()
    { return lowDimStencils_; }

    //! Return a reference to an empty stencil
    std::vector<unsigned int>& emptyStencil()
    { return emptyStencil_; }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    BulkProblem& bulkProblem_;
    LowDimProblem& lowDimProblem_;

    std::vector<BulkPointSource> bulkPointSources_;
    std::vector<LowDimPointSource> lowDimPointSources_;

    mutable std::vector<PointSourceData> pointSourceData_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimStencils_;
    std::vector<unsigned int> emptyStencil_;

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
    //! vector for the surface of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimSurfaceInBulkElement_;

    //! id generator for point sources
    unsigned int idCounter_;

    mutable std::map<std::pair<unsigned int, unsigned int>, BulkPrimaryVariables >  couplingSources_;
    mutable std::map<std::pair<unsigned int, unsigned int>, Scalar >  lowDimStencilLength_;
    mutable std::map<std::pair<unsigned int, unsigned int>, bool >  reuseCouplingSources_;
};

} // namespace Dumux

#endif
