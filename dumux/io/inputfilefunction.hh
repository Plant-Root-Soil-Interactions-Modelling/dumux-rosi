#ifndef INPUT_FILE_FUNCTION_HH
#define INPUT_FILE_FUNCTION_HH

#include <dumux/common/math.hh>
#include <assert.h>

namespace Dumux {

/**
 * Auxiliary class to pass a function via the input file, or via the grid file. The type is chosen automatically
 *
 * The type is determined from y, and x, that are read from the input file:
 *
 * Constant:    returns one constant value
 *              y is set to a single value, x is not set
 *
 * Table:       the value is calculated by linear interpolation from a table (y,x) given in the input file
 *              y and x are set to multiple sampling points
 *
 * Data:        the value is taken from the grid file
 *              y and x are not set
 *
 * Per Type     the value is given per type (e.g. soil layer, root order), where type is taken from the grid file.
 *              y is set to multiple values, x is not set (within the input file).
 *
 * Periodic     (optional). Periodic over one day, e.g. for root collar.
 *              y is set to two values determining min and max, x is not set
 *
 * If grid data is used (Data, or Per Type), the data must be set with InputFileFunction::setData, before calling InputFileFunction::f
 *
 */
class InputFileFunction {
public:

    enum {
        constant = 0, table = 1, data = 2, perType = 3, periodic = 4, empty = -1
    };

    //! don't know why it does not compile without it
    InputFileFunction() {
    }

    //!  constructor 1: handles one parameter of all types
    InputFileFunction(std::string nameY, std::string nameX, int dataIdx, int typeIdx = 0, bool enablePeriodic = false) {
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        try {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
                std::cout << "InputFileFunction: Constant (" << nameY << ")\n";
            } else {
                try {
                    xx_ = Dumux::getParam<std::vector<double>>(nameX);
                    type_ = table;
                    table_ = {xx_, yy_};
                    std::cout << "InputFileFunction: Table (" << nameY << ")\n";
                } catch (...) {
                    if ((enablePeriodic) && yy_.size() == 2) {
                        type_ = periodic;
                        std::cout << "InputFileFunction: Periodic (" << nameY << ")\n";
                    } else {
                        type_ = perType;
                        std::cout << "InputFileFunction: Constant per Type from Grid (" << nameY << ")\n";
                    }
                }
            }
        } catch (...) {
            type_ = data;
            std::cout << "InputFileFunction: Data from Grid (" << nameY << ")\n";
        }
    }

    //! constructor 2: class handles one parameter of type constant, perType, or data
    InputFileFunction(std::string nameY, int dataIdx, int typeIdx = 0) {
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        try {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
                std::cout << "InputFileFunction: Constant (" << nameY << ")\n";
            } else {
                type_ = perType;
                std::cout << "InputFileFunction: Constant per Type from Grid (" << nameY << ")\n";
            }
        } catch (...) {
            type_ = data;
            std::cout << "InputFileFunction: Data from Grid (" << nameY << ")\n";
        }
    }

    //! constructor 3: class handles one parameter of type constant, table, or periodic
    InputFileFunction(std::string nameY, std::string nameX, bool enablePeriodic = false) {
        yy_ = Dumux::getParam<std::vector<double>>(nameY);
        if (yy_.size() == 1) {
            type_ = constant;
            std::cout << "InputFileFunction: Constant (" << nameY << ")\n";
        } else {
            try {
                xx_ = Dumux::getParam<std::vector<double>>(nameX);
                type_ = table;
                table_ = {xx_, yy_};
                std::cout << "InputFileFunction: Table (" << nameY << ")\n";
            } catch (...) {
                if ((enablePeriodic) && yy_.size() == 2) {
                    type_ = periodic;
                    std::cout << "InputFileFunction: Periodic (" << nameY << ")\n";
                } else {
                    throw Dumux::ParameterException("InputFileFunction: Constructor 3, multiple values for parameter " + nameY);
                }
            }
        }
    }

    //! function
    double f(double x, size_t eIdx) const {
        switch (type_) {
        case constant: {
            return yy_[0];
        }
        case table: {
            return Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_);
        }
        case data: {
            return data_.at(eIdx);
        }
        case perType: {
            return yy_.at(size_t(data_.at(eIdx)));
        }
        case periodic: {
            return sin(x / (24. * 3600.) * 2. * M_PI) * (yy_[1] - yy_[0]) + 0.5 * (yy_[0] + yy_[1]);
        }
        default:
            throw Dumux::ParameterException("InputFileFunction: unknown function type");
        }
    }

    //! function (use for constructor 2)
    double f(size_t eIdx) const {
        assert((type_ != table && type_ != periodic && "InputFileFunction: call f(x) for table or periodic"));
        return f(0, eIdx);
    }

    //! function (use for constructor 3)
    double f(double x) const {
        assert((type_ != data && type_ != perType && "InputFileFunction: call f(eIdx) for data or perType"));
        return f(x, 0);
    }

    //! type
    int type() {
        return type_;
    }

    //! copies the data from a vector
    void setData(const std::vector<double> d) {
        if ((type_ == data) || (type_ == perType)) { // otherwise, don't bother
            data_ = d;
        }
    }

    // copies the data from grid file
    template<class GridData, class FVGridGeometry>
    void setGridData(const GridData& gridData, const FVGridGeometry& fvGridGeometry) {
        if ((type_ == data) || (type_ == perType)) { // otherwise, don't bother
            const auto& gridView = fvGridGeometry.gridView();
            data_.resize(gridView.size(0));
            const auto& elementMapper = fvGridGeometry.elementMapper();
            for (const auto& element : elements(gridView)) {
                const auto eIdx = elementMapper.index(element);
                const auto& elemParams = gridData.parameters(element);
                if (type_ == data) {
                    data_[eIdx] = elemParams[dataIdx_];
                } else { // type_ == perType
                    data_[eIdx] = elemParams[typeIdx_];
                }
            }
        }
    }

private:
    int type_ = empty;
    size_t dataIdx_ = -1;
    size_t typeIdx_ = -1;
    std::vector<double> xx_;
    std::vector<double> yy_;
    std::pair<std::vector<double>, std::vector<double>> table_;
    std::vector<double> data_;

};

} // namespace DUMUX

#endif
