#ifndef INPUT_FILE_FUNCTION_HH
#define INPUT_FILE_FUNCTION_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <assert.h>

namespace Dumux {

/**
 * Auxiliary class to pass a value or function in one argument via the input file, or via the grid file.
 * The type is chosen automatically.
 *
 * The type is determined from y, and x, that are read from the input file:
 *
 * Constant     returns one constant value
 *              y is set to a single value, x is not set
 *
 * Table        the value is calculated by linear interpolation from a table (y,x) given in the input file
 *              y and x are set to multiple sampling points
 *
 * Data         the value is taken from the grid file
 *              y and x are not set
 *
 * PerType      the value is constant and given per type (e.g. soil layer, root order), where type is taken from the grid file.
 *              y is set to multiple values, x is not set.

 * PerTypeIFF   the value is constant and given per type (e.g. soil layer, root order), where type is taken from another InputFileFunction calss
 *              y is set to multiple values, x is not set.
 *
 * TablePerType Multiple tables, named x0, x1, ... and y0, y1, ... The table number is chosen per type. (e.g. conductivities per root type)
 *              y and x are not set, x0,x1, ... , y0,y1, ... are set.
 *
 * Periodic     (optional). Periodic over one day, e.g. for root collar. // todo change this (all types can be periodic over the day)
 *              y is set to two values determining min and max, x is not set
 *
 * If grid data is used (Data, or Per Type), the data must be set with InputFileFunction::setData, before calling InputFileFunction::f
 */
class InputFileFunction {
public:

    enum {
        constant = 0, table = 1, data = 2, perType = 3, perTypeIFF = 4, periodic = 5, tablePerType = 6, empty = -1
    };

    //! don't know why it does not compile without it todo
    InputFileFunction() {
    }

    /**
     * Constructor 1: handles one parameter of all types
     *
     * todo passs std::string groupName, todo pass peridocity as parameter, todo move File for csv files here.
     *
     */
    InputFileFunction(std::string nameY, std::string nameX, int dataIdx, int typeIdx = 0, bool enablePeriodic = false, InputFileFunction* typeF = nullptr) {
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        nameY_ = nameY;
        bool hasX = Dumux::hasParam(nameY);
        bool hasY = Dumux::hasParam(nameY);
        if (hasY) {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
            } else {
                if (hasX) {
                    type_ = table;
                    xx_ = Dumux::getParam<std::vector<double>>(nameX);
                    table_.push_back( { xx_, yy_ });
                } else {
                    if ((enablePeriodic) && yy_.size() == 2) { // todo change periodicity
                        type_ = periodic;
                    } else {
                        if (typeF != nullptr) {
                            type_ = perTypeIFF;
                            iff_ = typeF;
                        } else {
                            type_ = perType;
                        }
                    }
                }
            }
        } else {
            try { // multiple tables todo remove try catch
                yy_ = Dumux::getParam<std::vector<double>>(nameY + std::to_string(0));
                table_.resize(0);
                int m = 0;
                while (true) {
                    try {
                        yy_ = Dumux::getParam<std::vector<double>>(nameY + std::to_string(m));
                        m++;
                    } catch (...) {
                        break;
                    }
                }
                for (int i = 0; i < m; i++) {
                    xx_ = Dumux::getParam<std::vector<double>>(nameX + std::to_string(i));
                    yy_ = Dumux::getParam<std::vector<double>>(nameY + std::to_string(i));
                    table_.push_back( { xx_, yy_ });
                }
                type_ = tablePerType;
            } catch (...) {
                type_ = data;
            }
        }
        cout();
    }

    //! constructor 2 (no tables): class handles one parameter of type constant, perType, or data
    InputFileFunction(std::string nameY, int dataIdx, int typeIdx = 0, InputFileFunction* typeF = nullptr) {
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        try { // todo remove try catch
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
            } else {
                if ((typeF != nullptr) && (typeF->type() == table)) { // todo: why only tables?
                    type_ = perTypeIFF;
                    iff_ = typeF;
                } else {
                    type_ = perType;
                }
            }
        } catch (...) {
            type_ = data;
        }
        cout();
    }

    //! constructor 3 (no grid data): class handles one parameter of type constant, table, or periodic
    InputFileFunction(std::string nameY, std::string nameX, bool enablePeriodic = false) {
        yy_ = Dumux::getParam<std::vector<double>>(nameY);
        if (yy_.size() == 1) {
            type_ = constant;
        } else {
            try { // todo remove try catch
                xx_ = Dumux::getParam<std::vector<double>>(nameX);
                type_ = table;
                table_.push_back( { xx_, yy_ });
            } catch (...) {
                if ((enablePeriodic) && yy_.size() == 2) {
                    type_ = periodic;
                } else {
                    throw Dumux::ParameterException("InputFileFunction: Constructor 3, multiple values for parameter " + nameY);
                }
            }
        }
        cout();
    }

    //! don't know why it does not compile without it
    InputFileFunction(std::vector<double> x, std::vector<double> y) {
        type_ = table;
        table_.resize(0);
        table_.push_back( { x, y });
        std::cout << "InputFileFunction: hard coded table \n";
    }


    /**
     * Evaluates the function
     *
     * @param x         the variable
     * @param eIdx      Dumux element index for grid data look up
     */
    double f(double x, size_t eIdx) const {
        switch (type_) {
        case constant: {
            return yy_[0];
        }
        case table: {
            return Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_[0]);
        }
        case data: {
            return data_.at(eIdx);
        }
        case perType: {
            return yy_.at(size_t(data_.at(eIdx)));
        }
        case perTypeIFF: {
            return yy_.at(size_t(iff_->f(x, eIdx) - 1)); // todo ?
        }
        case tablePerType: {
            size_t t = size_t(data_.at(eIdx)) - 1;
            assert( t>=0  && "type < 0" );
            assert( t<table_.size() && "type > read tables" );
            if(t>=table_.size() || (t<0) ) { // it seems assertions are not working ?????
                std::cout << "stranger things..." << t << std::flush;
            }
            return Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_.at(t));
        }
        case periodic: {
            double a = 0.5 * (yy_[1] - yy_[0]);
            return sin(x * 2. * M_PI - 0.5 * M_PI) * a + a + yy_[0];
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
    int type() const {
        return type_;
    }

    //! copies the data from a vector
    void setData(const std::vector<double> d) {
        if ((type_ == data) || (type_ == perType) || (type_ == tablePerType)) { // otherwise, don't bother
            data_ = d;
        }        std::cout << "InputFileFunction: hard coded table \n";
    }

    // copies the data from grid file
    template<class GridData, class FVGridGeometry>
    void setGridData(const GridData& gridData, const FVGridGeometry& fvGridGeometry) {
        if ((type_ == data) || (type_ == perType) || (type_ == tablePerType)) { // otherwise, don't bother
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

    /**
     * debugging output
     */
    void cout() {
        std::string s = "";
        if (periodic_) {
            s = " periodic over one day ";
        }
        switch (typeIdx_) {
        case(constant): std::cout << "InputFileFunction: Constant (" << nameY_ << ")"<< s << "\n"; break;
        case(table): std::cout << "InputFileFunction: Table (" << nameY_ << ")"<< s << "\n";  break;
        case(data): std::cout << "InputFileFunction: Data from Grid (" << nameY_ << ")"<< s << "\n"; break;
        case(perType): std::cout << "InputFileFunction: Constant per Type from Grid (" << nameY_ << ")"<< s << "\n"; break;
        case(perTypeIFF): std::cout << "InputFileFunction: Constant per Type from InputFileFunction (" << nameY_ << ")"<< s << "\n"; break;
        case(tablePerType): std::cout << "InputFileFunction: Table per Type from Grid (" << nameY_ << "), " << table_.size() << " types" << s << "\n"; break;
        default:
            std::cout << "InputFileFunction: unknown function type";
        }
    }

    int type_ = empty;
    size_t dataIdx_ = -1;
    size_t typeIdx_ = -1;
    std::string nameY_;
    bool periodic_ = false;
    InputFileFunction* iff_ = nullptr;
    std::vector<double> xx_;
    std::vector<double> yy_;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> table_;
    std::vector<double> data_;

};

} // namespace DUMUX

#endif
