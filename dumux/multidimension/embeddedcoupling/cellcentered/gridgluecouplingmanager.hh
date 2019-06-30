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
 *        domain. Intersection computation relies on the dune grid-glue backend.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_GRIDGLUE_EMBEDDEDCOUPLINGMANAGER_HH
#define DUMUX_CC_GRIDGLUE_EMBEDDEDCOUPLINGMANAGER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dumux/io/gridgluevtkwriter.hh>
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
 */
template<class TypeTag>
class CCGridGlueEmbeddedCouplingManager
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
    using BulkExtractor = typename Dune::GridGlue::Codim0Extractor<BulkGridView>;
    using LowDimExtractor = typename Dune::GridGlue::Codim0Extractor<LowDimGridView>;
    using GlueType = typename Dune::GridGlue::GridGlue<BulkExtractor, LowDimExtractor>;

public:

    /*!
     * \brief Constructor
     */
    CCGridGlueEmbeddedCouplingManager(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblem_(bulkProblem), lowDimProblem_(lowDimProblem)
    {
        // Check if we are using the cellcentered method in both domains
        static_assert(!bulkIsBox && !lowDimIsBox, "Using the cell-centered coupling manager for problems using box discretization!");
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
        asImp_().computePointSourceData(GET_RUNTIME_PARAM(TypeTag, int, Problem.IntegrationOrder),
                                        GET_RUNTIME_PARAM(TypeTag, bool, Problem.Verbose));
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
        // Initialize the descriptors (all elements)
        auto bulkDescriptor = [](const BulkElement& element, unsigned int subentity) { return true; };
        auto lowDimDescriptor = [](const LowDimElement& element, unsigned int subentity) { return true; };
        // construct element extractors
        auto bulkExtractor = std::make_shared<BulkExtractor>(this->bulkGridView(), bulkDescriptor);
        auto lowDimExtractor = std::make_shared<LowDimExtractor>(this->lowDimGridView(), lowDimDescriptor);
        // choose overlapping merger
        auto merger = std::make_shared<Dune::GridGlue::OverlappingMerge<bulkDim, lowDimDim, dimWorld>>();
        // build the glue object
        GlueType glue(bulkExtractor, lowDimExtractor, merger);

        // enable fallback brute-force step if the advancing front algorithm didn't find a seed
        // this is necessary of the grid has disconnected patches
        // unfortunately fallbacks for network grid will occur very often because the advancing
        // front algorithm in grid glue is not made for such grids. I might even occur
        // that the grid intersection fails for complicated networks. Hopefully there will be
        // a better implementation soon. See bounding box tree coupling manager for an alternative.
        merger->enableFallback(true);

        // compute intersections
        glue.build();

        // write output for debug purposes
        if (verbose)
        {
            GridGlueVtkWriter glueVtk;
            glueVtk.write(glue, "glue", verbose);
        }

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // a vector for the volumes the low-dim domain occupies in the bulk domain if it were full-dimensional
        lowDimVolumeInBulkElement_.resize(this->bulkGridView().size(0));

        // iterate over all intersections
        // we choose the target side (reverse, the low dimensional domain)
        for (const auto& is : intersections(glue, Dune::GridGlue::reversed))
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
            for (unsigned int outsideIdx = 0; outsideIdx < isIt->neighbor(0); ++outsideIdx)
            {
                const auto& outside = isIt->outside(outsideIdx);
                const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);
                if (lowDimDim == 1)
                    lowDimVolumeInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*M_PI*radius*radius;
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

    const BulkProblem& bulkProblem_;
    const LowDimProblem& lowDimProblem_;

    std::vector<BulkPointSource> bulkPointSources_;
    std::vector<LowDimPointSource> lowDimPointSources_;

    mutable std::vector<PointSourceData> pointSourceData_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimStencils_;
    std::vector<unsigned int> emptyStencil_;

    std::vector<Scalar> lowDimVolumeInBulkElement_;

    // id generator for point sources
    unsigned int idCounter_;
};

} // namespace Dumux

#endif
