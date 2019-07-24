#ifndef INPUT_FILE_FUNCTION_HH
#define INPUT_FILE_FUNCTION_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <assert.h>

#include "../../dumux/external/csv.h"

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
 *              y and x are both set to multiple sampling points
 *              y and x are not set, but a parameter File = example.csv, is set in the same parameter group.
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
 * All values options can be made sinusoidal over one day, by having the parameter sinusoidal = True in the same parameter group.
 *
 * If grid data is used (Data, or Per Type), the data must be set with InputFileFunction::setData, before calling InputFileFunction::f
 *
 */
class InputFileFunction {
public:

    enum {
        constant = 0, table = 1, data = 2, perType = 3, perTypeIFF = 4, tablePerType = 5, empty = -1
    };

    InputFileFunction() {  } // empty

    /**
     * Constructor 1: handles one parameter of all types
     *
     * call InputFileFunction::f(double x, size_t eIdx)
     */
    InputFileFunction(std::string groupName, std::string nameY, std::string nameX, int dataIdx, int typeIdx = 0, InputFileFunction* typeF = nullptr) {
        nameY = groupName +"." + nameY; // full names
        nameX = groupName +"." + nameX; // full names
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        nameY_ = nameY;
        sinusoidal_ = (hasParam(groupName +".Sinusoidal")) ? Dumux::getParam<bool>(groupName +".Sinusoidal") : false;
        if (Dumux::hasParam(nameY)) {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
            } else {
                if (Dumux::hasParam(nameX)) {
                    type_ = table;
                    xx_ = Dumux::getParam<std::vector<double>>(nameX);
                    table_.push_back( { xx_, yy_ });
                } else {
                    if (typeF != nullptr) {
                        type_ = perTypeIFF;
                        iff_ = typeF;
                    } else {
                        type_ = perType;
                    }
                }
            }
        } else {
            if (Dumux::hasParam(nameY + std::to_string(0))) {
                int m = 0;
                while (Dumux::hasParam(nameY + std::to_string(m))) { // count how many there are
                    m++;
                }
                for (int i = 0; i < m; i++) { // read them
                    xx_ = Dumux::getParam<std::vector<double>>(nameX + std::to_string(i));
                    yy_ = Dumux::getParam<std::vector<double>>(nameY + std::to_string(i));
                    table_.push_back( { xx_, yy_ });
                }
                type_ = tablePerType;
            } else {
                if  (hasParam(groupName +".File")) {
                    std::string name = Dumux::getParam<std::string>(groupName +".File");
                    std::string filename = getParam<std::string>(groupName +".File");
                    io::CSVReader<2> csv(filename);
                    csv.read_header(io::ignore_extra_column, "x", "y");
                    std::vector<double> x, y;
                    double a,b;
                    while(csv.read_row(a,b)){
                        x.push_back(a);
                        y.push_back(b);
                    }
                    table_.push_back( { x, y });
                    type_ = table;
                } else {
                    type_ = data;
                }
            }
        }
        cout();
    }

    /**
     * Constructor 2 (no tables):   handles one parameter of type constant, perType, perTypeIFF or data
     *                              i.e. no dependency of a variable
     *
     * a call to InputFileFunction::f(size_t eIdx) is sufficient
     */
    InputFileFunction(std::string groupName, std::string nameY, int dataIdx = 0, int typeIdx = 0, InputFileFunction* typeF = nullptr)
        :InputFileFunction(groupName, nameY, "no_valid_name", dataIdx, typeIdx, typeF) {
        assert(((type_==constant) || (type_==perType) || (type_==perTypeIFF) ||(type_==data)) &&
            "InputFileFunction: wrong type in constructor 2");
        cout();
    }

    /**
     * Constructor 3:   handles one parameter of all types but data, perType, TablePerType
     *                  i.e. all types without grid data
     *
     * a call to InputFileFunction::f(double x) is sufficient
     */
    InputFileFunction(std::string groupName, std::string nameY, std::string nameX, InputFileFunction* typeF = nullptr)
        :InputFileFunction(groupName, nameY, nameX, 0, 0, typeF) {
        assert(((type_==constant) || (type_==table) || (type_==perTypeIFF))  &&
            "InputFileFunction: wrong type in constructor 3");
        cout();
    }

    /*
     * Constructor 4:    a hard coded a look up table
     *
     * it should not be necessary to use this constructor
     * a call to InputFileFunction::f(double x) is sufficient
     */
    InputFileFunction(std::vector<double> x, std::vector<double> y) {
        type_ = table;
        table_.push_back( { x, y });
        std::cout << "InputFileFunction: hard coded table \n";
    }

    /**
     * Scales the dependent variable with @param s in InputFileFunction::f calls (for unit conversions)
     */
    void setVariableScale(double s) {
        vs_ = s;
    }

    /**
     * Scales the function's return value with @param s in InputFileFunction::f calls (for unit conversions)
     */
    void setFunctionScale(double s) {
        fs_ = s;
    }

    /**
     * Evaluates the function
     *
     * @param x         the variable (optionally is scaled, @see setVariableScale)
     * @param eIdx      Dumux element index for grid data look up
     *
     * @return          the function return values (optionally is scaled, @see setFunctionScale)
     */
    double f(double x, size_t eIdx) const {
        double fs  = fs_;
        if (sinusoidal_)  {
            fs = fs * (sin((x/86400.)*2.* M_PI - 0.5*M_PI) + 1);
            // x is assumed to be in seconds, scaling later, e.g. for look up tables based on days
        }
        x = x*vs_;

        switch (type_) {
        case constant: {
            return fs*yy_[0];
        }
        case table: {
            return fs*Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_[0]);
        }
        case data: {
            return fs*data_.at(eIdx);
        }
        case perType: {
            return fs*yy_.at(size_t(data_.at(eIdx)));
        }
        case perTypeIFF: {
            return fs*yy_.at(size_t(iff_->f(x, eIdx) - 1)); // we start at one, todo ?
        }
        case tablePerType: {
            size_t t = size_t(data_.at(eIdx)) - 1;
            assert( t>=0  && "InputFileFunction::f: table type < 0" );
            assert( t<table_.size() && "InputFileFunction::f: table type > available tables" );
            return Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_.at(t));
        }
        default:
            throw Dumux::ParameterException("InputFileFunction: unknown function type");
        }
    }

    /**
     *
     */
    double f(size_t eIdx) const {
        assert((type_ != table &&  "InputFileFunction: call f(x) for table"));
        return f(0, eIdx);
    }

    /**
     *
     */
    double f(double x) const {
        assert((type_ != data && type_ != perType && "InputFileFunction: call f(eIdx) for data or perType"));
        return f(x, 0);
    }

    /**
     * Type index of the InputFileFunction
     */
    int type() const {
        return type_;
    }

    /**
     * Copies the data from a vector, it is handled the same way as grid file data
     */
    void setData(const std::vector<double> d) {
        if ((type_ == data) || (type_ == perType) || (type_ == tablePerType)) { // otherwise, don't bother
            data_ = d;
        }
    }

    /**
     * Copies the data from the grid file
     */
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

    /**
     * Quick info about the type of function
     */
    void cout() {
        std::string s = "";
        if (sinusoidal_) {
            s = " sinusoidal over one day ";
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

private:

    int type_ = empty;
    size_t dataIdx_ = -1;
    size_t typeIdx_ = -1;
    std::string nameY_;
    bool sinusoidal_ = false;
    double vs_ = 1.; // for unit conversions
    double fs_ = 1.; // for unit conversions
    InputFileFunction* iff_ = nullptr;
    std::vector<double> xx_;
    std::vector<double> yy_;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> table_;
    std::vector<double> data_;

};

} // namespace DUMUX

#endif
